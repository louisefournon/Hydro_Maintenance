#include <LagBFunction.h>
#include <BendersBFunction.h>
#include <BendersBlock.h>
#include <AbstractBlock.h>
#include <Block.h>
#include <ColVariable.h>
#include <FRowConstraint.h>
#include <FRealObjective.h>
#include <UCBlock.h>
#include <IntermittentUnitBlock.h>
#include <HydroUnitMaintenance.h>
#include <HydroUnitBlock.h>
#include <Function.h>
#include <Block.h>


using namespace SMSpp_di_unipi_it;

std::string filename{};

int main(){
    
    ////////////////////////////////////////////////////////////////////////////////
    //                             Data extraction                                //
    ////////////////////////////////////////////////////////////////////////////////
    
    
    // Read the netCDF File with all the data

    netCDF::NcFile f;
    try {
        f.open( "20_10_h_w.nc4", netCDF::NcFile::read );
    } catch( netCDF::exceptions::NcException & e ) {
        std::cerr << "Cannot open nc4 file " << filename << std::endl;
        exit( 1 );
    }

    netCDF::NcGroupAtt gtype = f.getAtt( "SMS++_file_type" );
    if( gtype.isNull() ) {
        std::cerr << filename << " is not an SMS++ nc4 file" << std::endl;
        exit( 1 );
    }

    int type;
    gtype.getValues( &type );

    if( type != eBlockFile ) {
        std::cerr << filename << " is not an SMS++ nc4 Block file" << std::endl;
        exit( 1 );
    }
    
    // Get the data of the general UC block (bg) contained in Block_0

    netCDF::NcGroup bg = f.getGroup( "Block_0" );

    if( bg.isNull() ) {
        std::cerr << "Block_0 empty or undefined in " << filename << std::endl;
        exit( 1 );
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    //                 Construction of the initial problem W(z = 0) as defined in section 4                //
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // In the SMS++ modelization, W(z) is composed of a "general Block" representing the overall problem and constraints (demand constraints etc.) and several thermal and hydro blocks representing our assets.
    // Note that in this first modelization, W(z) doesn't actually depends on z : the variable xi created by the splitting of the maintenance variable z is set as equal to 0 (i.e. xi = 0 and not xi = z). That constraint will be updated to xi = z later in the model, in the BendersBFunction. 
    
    // 1) Construct the UCBlock representing W(z)

    auto uc_block = dynamic_cast<UCBlock *>(Block::new_Block( "UCBlock" ));
    uc_block->deserialize( bg );
    
    
    // 1.0) Get dimensions of our constraints
    
    auto network_data = uc_block->get_NetworkData();
    unsigned int number_nodes = network_data ? network_data->get_number_nodes() : 1;
    unsigned int time_horizon = uc_block->get_time_horizon();
    unsigned int nb_primary_zones = uc_block->get_number_primary_zones();
    unsigned int nb_secondary_zones = uc_block->get_number_secondary_zones();

    
    ////////////////////////////////////////////////////////////////////////////////
    //       Construction of the lagrangian of W(0) : constraints relaxation      //
    ////////////////////////////////////////////////////////////////////////////////
    
    
    // 1.1) Relax the desired constraints in uc_block.
    
    
    // First we relax the demand constraints which are the first 3 constraints: v_node_injection_constraints, v_PrimaryDemand_Const, v_SecondaryDemand_Const 

    auto node_injection_constraints = uc_block->get_static_constraint<FRowConstraint, 2>(1);  
    auto primary_demand_constraints = uc_block->get_static_constraint<FRowConstraint, 2>(2); 
    auto secondary_demand_constraints = uc_block->get_static_constraint<FRowConstraint, 2>(3); 
    
//    un_any_static( node_injection_constraints , []( Constraint * constraint ) { constraint->relax( true ); } , un_any_type<Constraint> );
//    un_any_static( primary_demand_constraints , []( Constraint * constraint ) { constraint->relax( true ); } , un_any_type<Constraint> );
//    un_any_static( secondary_demand_constraints , []( Constraint * constraint ) { constraint->relax( true ); } , un_any_type<Constraint> );
    
    // Then we look for the constraint x = \xi in each nested block to relax them
    // We also want the indexes and number of hydroUnitBlocks for the rest
    
    std::vector<int> idx_hydro_blocks; 
    int nb_hydro_blocks;
    
    auto sb =  uc_block->get_nested_Blocks();
    
    for( UnitBlock::Index i = 0 ; i < sb.size() ; ++i ) { // Each i is a subblock (thermal or hydro)
        
        auto unit_block = dynamic_cast<UnitBlock *>( sb[i] );
        auto hydro_unit_block = dynamic_cast<HydroUnitBlock *>(unit_block); // I believe this only creates a block if it's a hydroUnitBlock right ?
        if( hydro_unit_block != nullptr ){
            idx_hydro_blocks.push_back(i);
            auto constraints = hydro_unit_block->get_static_constraint<FRowConstraint, 2>(2);  // "Splitting var" is the second constraint
//            un_any_static( constraints , []( Constraint * constraint ) { constraint->relax( true ); } , un_any_type<Constraint> );
        }
    }
    
    nb_hydro_blocks = idx_hydro_blocks.size();
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //          Construction of the object containing the whole lagrangian problem (i.e. lagrangian_block)          //
    //                                 whose supremum will be the theta_0 function                                  //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Note: We have 2 distincts objects here: lagrangian_block and lagrangian_function. lagrangian_block is the whole lagrangian problem, and it contains lagrangian_function (it's its objective function). The dual function will be obtained after solving the lagrangian_block.
    
    // 2) Construct the Lagrangian Block (the Block whose objective will be a LagBFunction).

    auto lagrangian_block = new AbstractBlock();
    
    
    ///////////////////////////////////////////////////////////////////////////////////////
    //              Creation of lambda (dual variables) for the lagrangian               //
    ///////////////////////////////////////////////////////////////////////////////////////
    

        //2.1) Create the variables to the lagrangian_block:
    
    // We have 3 lambdas for the 3 demand constraints (dimension = 1) and nb_hydro_blocks other lambdas for the splittingVar constraints (dimension = nb_generators = f_number_arcs)
    
    // Create the dual variables for the demand constraints
    
    std::vector< std::vector< ColVariable > > lambda_node_injection(time_horizon, std::vector< ColVariable >(number_nodes));
    std::vector< std::vector< ColVariable > > lambda_primary_demand(time_horizon, std::vector< ColVariable >(nb_primary_zones));
    std::vector< std::vector< ColVariable > > lambda_secondary_demand(time_horizon, std::vector< ColVariable >(nb_secondary_zones));
    
    // Create the lambdas of the Splitting variables constraints (xi = z)
    
    std::vector<UnitBlock::Index> nb_generators(nb_hydro_blocks, 0.0); // Number of generators for each hydroUnit (assets)
    
    std::vector<std::vector<std::vector< ColVariable > > > lambda_splitting_var; 
    
    for( UnitBlock::Index i : idx_hydro_blocks ) {
    
        auto unit_block = dynamic_cast<UnitBlock *>( sb[i] );
        auto hydro_unit_block = dynamic_cast<HydroUnitBlock *>(unit_block);
        
        if( hydro_unit_block != nullptr ){
            UnitBlock::Index n_i = hydro_unit_block->get_number_generators();
            std::vector< std::vector<ColVariable>> subLambda(time_horizon, std::vector< ColVariable >(n_i));
            lambda_splitting_var.push_back(subLambda); 
            nb_generators[i] = n_i;
        }
    }
    
        //2.2) Set the sign of the variables if necessary:
    
    for( int t = 0; t < time_horizon; t++){
        for( int j = 0; j < number_nodes; j++){
            lambda_node_injection[t][j].is_positive( true );
        }
        for( int j = 0; j < nb_primary_zones; j++){
            lambda_primary_demand[t][j].is_positive( true );
        }
        for( int j = 0; j < nb_secondary_zones; j++){
            lambda_secondary_demand[t][j].is_positive( true );
        }
    }
    
    for(int i = 0; i < nb_hydro_blocks; i++) {
        for(int t = 0; t < time_horizon; t++){
            for(UnitBlock::Index j = 0; j < nb_generators[i]; j++){
                lambda_splitting_var[i][t][j].is_positive( true );
            }
        }
    }

        //2.3) Add the variables to the Lagrangian Block:

    lagrangian_block->add_static_variable( lambda_node_injection );
    lagrangian_block->add_static_variable( lambda_primary_demand );
    lagrangian_block->add_static_variable( lambda_secondary_demand );
    lagrangian_block->add_static_variable( lambda_splitting_var );
    
    ///////////////////////////////////////////////////////////////////////////////////
    //       Construction of the actual Lagrangian function (lagrangian_function)    //
    ///////////////////////////////////////////////////////////////////////////////////
    

        //2.4) Construct the LagBFunction (objective of our lagrangian_block):

    LagBFunction lagrangian_function;
    lagrangian_function.set_inner_block( uc_block );
    
    
    

        //2.5) Associate each lambda_i with a relaxed function (a Function belonging to a relaxed RowConstraint; FRowConstraint has the method get_function() to retrieve a pointer to the Function associated with that Constraint):
    
    for( int t = 0; t < time_horizon; t++){

        for( int j = 0; j < number_nodes; j++){
            
            Function * function = node_injection_constraints[t][j].get_function();
            std::pair < std::vector< ColVariable >, Function * > node_injection( lambda_node_injection[t][j], function );
            lagrangian_function.set_dual_pairs(node_injection);
        }
        for( int j = 0; j < nb_primary_zones; j++){
            
            Function * function = primary_demand_constraints[t][j].get_function(); 
            std::pair < std::vector< ColVariable >, Function * > primary_demand( lambda_primary_demand[t][j],function );

            lagrangian_function.set_dual_pairs(primary_demand);
        }
        for( int j = 0; j < nb_secondary_zones; j++){
            
            Function * function = secondary_demand_constraints[t][j].get_function();
            std::pair < std::vector< ColVariable >, Function * > secondary_demand( lambda_secondary_demand[t][j], function );
            lagrangian_function.set_dual_pairs(secondary_demand);
        }
    }
    
    
    for(int i = 0; i < nb_hydro_blocks; i++) {
        
        for(int t = 0; t < time_horizon; t++){
            
            for(UnitBlock::Index j = 0; j < nb_generators[i]; j++){
                
                Function * function = splitting_var_constraints[t][j].get_function();
                std::pair < std::vector< ColVariable >, Function * > splitting_var( lambda_splitting_var[i][t][j], function );
                lagrangian_function.set_dual_pairs(splitting_var);
                
            }
        }
    }

        //2.6) Add the Objective to the Lagrangian Block.

    auto objective = new FRealObjective( lagrangian_block , & lagrangian_function );
    objective->set_sense( Objective::eMax );
    lagrangian_block->set_objective( objective );
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //            Creation of the \bar W(z) function = supremum of theta_0 updated to theta_z using BendersBlock             //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    

    //3) If you want, construct a BendersBlock that will have the BendersBFunction as the objective function. Let m be the number of z variables i.e. the total number of generators in our system 
    

    BendersBlock benders_block( nullptr , m);

    //4) Construct the BendersBFunction.

        //4.1) Create a vector of pointers to z (z doesn't exist yet in our UCBlock, for now we have \xi = 0. BendersBFunction will later update the constraints of UCBlock using the mapping to include the dependance in z in the constraints \xi = z):

    std::vector< ColVariable * > z( m );
    for( auto & z_i : benders_block.get_variables() )
        z.push_back( const_cast<ColVariable *>( & z_i ) );
    
    // Creation of the mapping that will allow us to update our theta_0 problem to theta_z
    
        //4.2) Create the mapping formed by the m x m identity matrix A, the m-dimensional zero vector b, the m-dimensional ConstraintSide vector with all components equal to BendersBFunction::ConstraintSide::eBoth, and the vector of pointers to the RowConstraints representing z_i = \xi_i that you defined in the HydroUnitBlock:
    
    // Create identity matrix
    std::vector< std::vector<Function::FunctionValue> > A(m, std::vector<Function::FunctionValue>(m,0));
    for(unsigned int t = 0; t < m; t++)
        A[t][t] = 1;

    // Create b vector
    std::vector<Function::FunctionValue> b(m,0);
    
    // Create sides vector
    std::vector<BendersBFunction::ConstraintSide> sides(m, BendersBFunction::ConstraintSide::eBoth);
    
    // Create constraints pointers vector
    std::vector< RowConstraint * > constraints_pointer;
    
    for( UnitBlock::Index i : idx_hydro_blocks ) {
        
        auto unit_block = dynamic_cast<UnitBlock *>( sb[i] );
        auto hydro_unit_block = dynamic_cast<HydroUnitBlock *>(unit_block);
        if( hydro_unit_block != nullptr ){
//            auto constraint = hydro_unit_block->get_static_constraint(1);
//            constraints_pointer.push_back( &constraint);
        }
    }

    BendersBFunction benders_function( lagrangian_block ,  std::move(z) ,  std::move(A) ,  std::move(b) , std::move(constraints_pointer) , std::move(sides) );

    //5) Add the BendersBFunction as the Objective of the BendersBlock

    benders_block.set_function( & benders_function );
    
    
}