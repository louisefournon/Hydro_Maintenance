#include <LagBFunction.h>
#include <BendersBFunction.h>
#include <BendersBlock.h>
#include <AbstractBlock.h>
#include <Block.h>
#include <ColVariable.h>
#include <FRowConstraint.h>
#include <FRealObjective.h>

//// Test GIT

using namespace SMSpp_di_unipi_it;


int main(){
    
    // Read the netCDF File with all the data

    netCDF::NcFile file_path;
    try {
        f.open( filename, netCDF::NcFile::read );
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
    
    // Get the data of the general UC block (bg)

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

    auto linking_const = uc_block->get_static_constraints();    
    linking_const(v_node_injection_constraints).relax(true);
    linking_const(v_PrimaryDemand_Const).relax(true);
    linking_const(v_SecondaryDemand_Const).relax(true); // I'm not sure this is the right way to access the constraints
    
    // Then we look for the constraint x = \xi in each nested block to relax them
    
    // We also want the indexes and number of hydroUnitBlocks for the rest
    vector<int> idx_hydro_blocks; 
    int nb_hydro_blocks;
    
    for( auto i: uc_block->get_nested_Blocks() ) {
        
        auto unit_block = dynamic_cast<UnitBlock *>(i);
        auto hydro_unit_block = dynamic_cast<HydroUnitBlock *>(unit_block); // I believe this only creates a block if it's a hydroUnitBlock right ?
        if( hydro_unit_block != nullptr ){
            idx_hydro_blocks.push_back(i);
            auto constraint = hydro_unit_block->get_static_constraints("XiEqualZ");
            constraint.relax(true);
        }
    }
    
    nb_hydro_blocks = idx_hydro_blocks.size();
    
    // 2) Construct the Lagrangian Block (the Block whose objective will be a LagBFunction).

    auto lagrangian_block = new AbstractBlock();

        //2.1) Create the variables to the lagrangian_block:
    
    // We have 3 lambdas for the 3 demand constraints (dimension = 1) and nb_hydro_blocks other lambdas for the XiEqualZ constraints (dimension = nb_generators)
    
    // Create the dual variables for the demand constraints
    std::vector< ColVariable > lambda_1(1);
    std::vector< ColVariable > lambda_2(1);
    std::vector< ColVariable > lambda_3(1);
    
    // Create the lambdas of the XiEqualZ constraints 
    
    std::vector< std::vector< ColVariable > > lambda;
    for( auto i: uc_block->get_nested_Blocks() ) {
        
        auto unit_block = dynamic_cast<UnitBlock *>(i);
        auto hydro_unit_block = dynamic_cast<HydroUnitBlock *>(unit_block);
        if( hydro_unit_block != nullptr ){
            lambda.push_back(std::vector< ColVariable > subLambda(hydro_unit_block.f_number_arcs)); // Not sure about that attribute access
            auto constraint = hydro_unit_block->get_static_constraints("XiEqualZ");
            constraint.relax(true);
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // I haven't done the rest yet
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //2.2) Set the sign of the variables if necessary:

    for( auto & lambda_i : lambda ) {
        lambda_i.is_positive( true );
    }

        //2.3) Add the variables to the Lagrangian Block:

    lagrangian_block->add_static_variable( lambda );

        //2.4) Construct the LagBFunction:

    LagBFunction lagrangian_function;
    lagrangian_function.set_inner_block( & uc_block );

        //2.5) Associate each lambda_i with a relaxed function (a Function belonging to a relaxed RowConstraint; FRowConstraint has the method get_function() to retrieve a pointer to the Function associated with that Constraint):

    lagrangian_function.set_dual_pairs( ... );

        //2.6) Add the Objective to the Lagrangian Block.

    auto objective = new FRealObjective( lagrangian_block , & lagrangian_function );
    objective->set_sense( Objective::eMax );
    lagrangian_block->set_objective( objective );

    //3) If you want, construct a BendersBlock that will have the BendersBFunction as the objective function. Let n be the number of z variables.

    BendersBlock benders_block( nullptr , n );

    //4) Construct the BendersBFunction.

        //4.1) Create a vector of pointers to z:

    std::vector< ColVariable * > z( n );
    for( auto & z_i : benders_block.get_variables() )
        z.push_back( const_cast<ColVariable *>( & z_i ) );

        //4.2) Create the mapping formed by the n x n identity matrix A, the n-dimensional zero vector b, the n-dimensional ConstraintSide vector with all components equal to BendersBFunction::ConstraintSide::eBoth, and the vector of pointers to the RowConstraints representing z_i = \xi_i that you defined in the HydroUnitBlock:

    BendersBFunction benders_function( lagrangian_block , z , A , b , constraints , sides );

    //5) Add the BendersBFunction as the Objective of the BendersBlock

    benders_block.set_function( & benders_function );
    
    
}