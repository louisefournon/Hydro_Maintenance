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


using namespace SMSpp_di_unipi_it;

std::string filename{};

int main(){
    
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

    
    // 1) Construct the UCBlock.

    auto uc_block = dynamic_cast<UCBlock *>(Block::new_Block( "UCBlock" ));
    uc_block->deserialize( bg );
    
    // 1.1) Relax the desired constraints in uc_block.
        // Constraint has a method called relax(). You should invoke this method for each Constraint you want to relax (and pass true as argument).
    
    
    // First we relax the demand constraints which are the first 3 constraints: v_node_injection_constraints, v_PrimaryDemand_Const, v_SecondaryDemand_Const 

//    auto linking_const_0 = uc_block->get_static_constraint(0);  
//    auto linking_const_1 = uc_block->get_static_constraint(1); 
//    auto linking_const_2 = uc_block->get_static_constraint(2); 
//    
//    un_any_static( linking_const_0 , []( Constraint * constraint ) { constraint->relax( true ); } , un_any_type<Constraint> );
//    un_any_static( linking_const_1 , []( Constraint * constraint ) { constraint->relax( true ); } , un_any_type<Constraint> );
//    un_any_static( linking_const_2 , []( Constraint * constraint ) { constraint->relax( true ); } , un_any_type<Constraint> );
    
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
//            auto constraints = hydro_unit_block->get_static_constraint(1);  // "XiEqualZ" is the second constraint
//            un_any_static( constraints , []( Constraint * constraint ) { constraint->relax( true ); } , un_any_type<Constraint> );
        }
    }
    
    nb_hydro_blocks = idx_hydro_blocks.size();
    
    // 2) Construct the Lagrangian Block (the Block whose objective will be a LagBFunction).

    auto lagrangian_block = new AbstractBlock();

        //2.1) Create the variables to the lagrangian_block:
    
    // We have 3 lambdas for the 3 demand constraints (dimension = 1) and nb_hydro_blocks other lambdas for the XiEqualZ constraints (dimension = nb_generators = f_number_arcs)
    
    // Create the dual variables for the demand constraints
    std::vector< ColVariable > lambda_0(1);
    std::vector< ColVariable > lambda_1(1);
    std::vector< ColVariable > lambda_2(1);
    
    // Create the lambdas of the XiEqualZ constraints 
    
    int m = 0; // Total number of generators (assets)
    
    std::vector< std::vector< ColVariable > > lambda;
    for( UnitBlock::Index i : idx_hydro_blocks ) {
        
        auto unit_block = dynamic_cast<UnitBlock *>( sb[i] );
        auto hydro_unit_block = dynamic_cast<HydroUnitBlock *>(unit_block);
        if( hydro_unit_block != nullptr ){
            UnitBlock::Index n_i = hydro_unit_block->get_number_generators();
            std::vector< ColVariable > subLambda(n_i);
            lambda.push_back(subLambda); 
            m += n_i;
        }
    }
    
    
        //2.2) Set the sign of the variables if necessary:
    
    for( int i = 0; i < lambda_0.size(); i++){
        lambda_0[i].is_positive( true );
    }
    for( int i = 0; i < lambda_1.size(); i++){
        lambda_1[i].is_positive( true );
    }
    for( int i = 0; i < lambda_2.size(); i++){
        lambda_2[i].is_positive( true );
    }

    for( auto & lambda_i : lambda ) {
        for(int k = 0; k < lambda_i.size(); k++){
            lambda_i[k].is_positive( true );
        }
    }

        //2.3) Add the variables to the Lagrangian Block:

    lagrangian_block->add_static_variable( lambda );
    lagrangian_block->add_static_variable( lambda_0 );
    lagrangian_block->add_static_variable( lambda_1 );
    lagrangian_block->add_static_variable( lambda_2 );


        //2.4) Construct the LagBFunction:

    LagBFunction lagrangian_function;
    lagrangian_function.set_inner_block( uc_block );
    
    

        //2.5) Associate each lambda_i with a relaxed function (a Function belonging to a relaxed RowConstraint; FRowConstraint has the method get_function() to retrieve a pointer to the Function associated with that Constraint):
    
//    std::pair l_0( lambda_0, linking_const[0].get_function() );
//    std::pair l_1( lambda_1, linking_const[1].get_function() );
//    std::pair l_2( lambda_2, linking_const[2].get_function() );
//                  
//    
//    lagrangian_function.set_dual_pairs( l_0 );
//    lagrangian_function.set_dual_pairs( l_1 );
//    lagrangian_function.set_dual_pairs( l_2 );
    
    for( UnitBlock::Index i : idx_hydro_blocks ) {
        
        auto unit_block = dynamic_cast<UnitBlock *>( sb[i] );
        auto hydro_unit_block = dynamic_cast<HydroUnitBlock *>(unit_block);
        if( hydro_unit_block != nullptr ){
//            auto constraint = hydro_unit_block->get_static_constraint(1);
//            lagrangian_function.set_dual_pairs( std::pair l( lambda[i], constraint.get_function() ) );
        }
    }

        //2.6) Add the Objective to the Lagrangian Block.

    auto objective = new FRealObjective( lagrangian_block , & lagrangian_function );
    objective->set_sense( Objective::eMax );
    lagrangian_block->set_objective( objective );

    //3) If you want, construct a BendersBlock that will have the BendersBFunction as the objective function. Let m be the number of z variables i.e. the total number of generators in our system 
    

    BendersBlock benders_block( nullptr , m);

    //4) Construct the BendersBFunction.

        //4.1) Create a vector of pointers to z (z doesn't exist yet in our UCBlock, for now we have \xi = 0. BendersBFunction will later update the constraints of UCBlock using the mapping to include the dependance in z in the constraints \xi = z):

    std::vector< ColVariable * > z( m );
    for( auto & z_i : benders_block.get_variables() )
        z.push_back( const_cast<ColVariable *>( & z_i ) );
    
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
    
    
    

    
    
    

    BendersBFunction benders_function( lagrangian_block , z , A , b , constraints_pointer , sides );

    //5) Add the BendersBFunction as the Objective of the BendersBlock

    benders_block.set_function( & benders_function );
    
    
}