#include <LagBFunction.h>
#include <BendersBFunction.h>
#include <BendersBlock.h>
#include <AbstractBlock.h>
#include <Block.h>
#include <ColVariable.h>
#include <FRowConstraint.h>
#include <FRealObjective.h>

using namespace SMSpp_di_unipi_it;


int main(){
    
    // 1) Construct the UCBlock.

    auto uc_block = Block::deserialize( file_path );
    
    for( auto i: uc_block->get_nested_Blocks() ) {
        
        auto unit_block = dynamic_cast<UnitBlock *>(i);
        auto hydro_unit_block = dynamic_cast<HydroUnitBlock *>(unit_block);

        auto constraint = hydro_unit_block->get_static_constraints();

                auto active_power = hydro_unit_block->get_active_power( g );
                std::cout << "active_power     = [";
                for( UnitBlock::Index t = 0; t < unit_block->get_time_horizon(); ++t ) {
                    std::cout << std::setw( 20 ) << active_power[t].get_value();
                    }
                std::cout << " ]" << std::endl;
                }

        }

    }

        // 1.1) Relax the desired constraints in uc_block.
        // Constraint has a method called relax(). You should invoke this method for each Constraint you want to relax (and pass true as argument).

    // 2) Construct the Lagrangian Block (the Block whose objective will be a LagBFunction).

    auto lagrangian_block = new AbstractBlock();

        //2.1) Create the variables to the lagrangian_block:

    std::vector< ColVariable > lambda( size_of_lambda );

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