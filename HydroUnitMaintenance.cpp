/*--------------------------------------------------------------------------*/
/*--------------------- File HydroUnitMaintenance.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the HydroUnitMaintenance class.
 *
 * \version 0.11
 *
 * \date 30 - 03 - 2020
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 *
 * \author Ali Ghezelsoflu \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 *
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <random>
#include "HydroUnitMaintenance.h"
#include "LinearFunction.h"
#include <map>
#include "FRowConstraint.h"
#include "UnitBlock.h"


/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;


/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register HydroUnitMaintenance to the Block factory


SMSpp_insert_in_factory_cpp_1( HydroUnitMaintenance );

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS OF HydroUnitMaintenance -----------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
void HydroUnitMaintenance::deserialize( netCDF::NcGroup & group ) {

 UnitBlock::deserialize_time_horizon( group );
 UnitBlock::deserialize_change_intervals( group );

 ::deserialize_dim( group, "NumberReservoirs", f_number_reservoirs, true );
 ::deserialize_dim( group, "NumberArcs", f_number_arcs, true );
 f_number_reservoirs = f_number_reservoirs ? f_number_reservoirs : 1;
 f_number_arcs = f_number_arcs ? f_number_arcs : 1;

 ::deserialize( group, "NumberPieces", f_number_arcs,
                v_number_pieces, true, true );


 for( auto & n : v_number_pieces ) {
  f_total_number_pieces += n;
 }

 f_total_number_pieces = f_total_number_pieces ?
                         f_total_number_pieces : f_number_arcs;

 ::deserialize( group, "StartArc", f_number_arcs, v_start_arc );
 ::deserialize( group, "EndArc", f_number_arcs, v_end_arc );

 ::deserialize( group, "Inflows", v_inflows, true, false );
 transpose( v_inflows );

 ::deserialize( group, "MinFlow", v_minimum_flow, true, true );
 transpose( v_minimum_flow );

 ::deserialize( group, "MaxFlow", v_maximum_flow, true, true );
 transpose( v_maximum_flow );


 ::deserialize( group, "MinPower", v_minimum_power, true, true );
 transpose( v_minimum_power );

 ::deserialize( group, "MaxPower", v_maximum_power, true, true );
 transpose( v_maximum_power );

 ::deserialize( group, "DeltaRampUp", v_delta_ramp_up, true, true );
 transpose( v_delta_ramp_up );

 ::deserialize( group, "DeltaRampDown", v_delta_ramp_down, true, true );
 transpose( v_delta_ramp_down );

 ::deserialize( group, "PrimaryRho", v_primary_rho, true, true );
 transpose( v_primary_rho );

 ::deserialize( group, "SecondaryRho", v_secondary_rho, true, true );
 transpose( v_secondary_rho );

 ::deserialize( group, "LinearTerm", f_total_number_pieces,
                v_linear_term, true, true );

 ::deserialize( group, "ConstantTerm", f_total_number_pieces,
                v_const_term, true, true );

 ::deserialize( group, "InertiaPower", v_inertia_power, true, true );

 ::deserialize( group, "InitialFlowRate", f_number_arcs,
                v_initial_flow_rate, true, true );

 ::deserialize( group, "InitialVolumetric", f_number_reservoirs,
                v_initial_volumetric, true, true );

 ::deserialize( group, "UphillFlow", f_number_arcs,
                v_uphill_delay, true, true );

 ::deserialize( group, "DownhillFlow", f_number_arcs,
                v_downhill_delay, true, true );

 ::deserialize( group, "MinVolumetric", v_minimum_volumetric, true, true );
 transpose( v_minimum_volumetric );
 ::deserialize( group, "MaxVolumetric", v_maximum_volumetric, true, true );
 transpose( v_maximum_volumetric );

 decompress_array( v_minimum_flow );
 decompress_array( v_maximum_flow );
 decompress_vol( v_minimum_volumetric );
 decompress_vol( v_maximum_volumetric );
 //decompress_vol( v_inflows );

 decompress_array( v_minimum_power );
 decompress_array( v_maximum_power );
 decompress_array( v_delta_ramp_up );
 decompress_array( v_delta_ramp_down );
 decompress_array( v_primary_rho );
 decompress_array( v_secondary_rho );
 decompress_array( v_inertia_power );

 if( v_linear_term.size() == 1 ) {
  v_linear_term.resize( f_number_arcs, v_linear_term[ 0 ] );
 }
 if( v_const_term.size() == 1 ) {
  v_const_term.resize( f_number_arcs, v_const_term[ 0 ] );
 }

 UnitBlock::deserialize( group );
}// end( HydroUnitMaintenance::deserialize )

/*--------------------------------------------------------------------------*/

void HydroUnitMaintenance::generate_abstract_variables( Configuration *stvv )
{
 UnitBlock::generate_abstract_variables( stvv );

 if( f_time_horizon == 0 ) {
  // there are no variables to be generated
  return;
 }

 if( !v_volumetric.empty()  ) {
  // the abstract variables should be generated only once
  return;
 }
 v_volumetric.resize(boost::extents[f_number_reservoirs][ f_time_horizon]);
 
 add_static_variable ( v_volumetric, "vol" );

 if( !v_flow_rate.empty() ) {
  // the abstract variables should be generated only once
  return;
 }
 v_flow_rate.resize(boost::extents[f_time_horizon][f_number_arcs]);
 for( Index t = 0; t < f_time_horizon; ++t ) {
  for( Index g = 0; g < f_number_arcs; ++g ) {
   auto & flow_rate = v_flow_rate[ t ][ g ];

   flow_rate.set_type( ColVariable::kContinuous );
  }
 }
 add_static_variable ( v_flow_rate, "F" );

 if( !v_active_power.empty()  ) {
  // the abstract variables should be generated only once
  return;
 }

 v_active_power.resize(boost::extents[f_time_horizon][f_number_arcs]);

 for( Index t = 0; t < f_time_horizon; ++t ) {
 for( Index g = 0; g < f_number_arcs; ++g ) {
   v_active_power[ t ][ g ].set_type( ColVariable::kContinuous );

  }
 }
 add_static_variable ( v_active_power, "p" );

    
/////////////////////////////////////////////////   
//       Add the maintenance variables         //
/////////////////////////////////////////////////
    
 if( !v_maintenance.empty()  ) {
  // the abstract variables should be generated only once
  return;
 }

 v_maintenance.resize(boost::extents[f_time_horizon][f_number_arcs]);

 for( Index t = 0; t < f_time_horizon; ++t ) {
     for( Index g = 0; g < f_number_arcs; ++g ) {
        v_maintenance[ t ][ g ].set_type( ColVariable::kContinuous ); // v_maintenance = 0 if maintenance, 1 otherwise
        // Note that since we have v_maintenance = z, we don't have to specify that v_maintenance is binary

        }
 }
 add_static_variable ( v_maintenance, "z" );
    
//////////////////////////////////////////////////////////
    
    

  if( !v_primary_spinning_reserve.empty()  ) {
   // the abstract variables should be generated only once
   return;
  }
  v_primary_spinning_reserve.resize(boost::extents[f_time_horizon][f_number_arcs]);
 for( Index t = 0; t < f_time_horizon; ++t ) {
 for( Index g = 0; g < f_number_arcs; ++g ) {
    v_primary_spinning_reserve[ t ][ g ].set_type( ColVariable::kNonNegative );

   }
  }
  add_static_variable ( v_primary_spinning_reserve, "pr" );

  if( !v_secondary_spinning_reserve.empty()  ) {
   // the abstract variables should be generated only once
  return;
  }
 v_secondary_spinning_reserve.resize(boost::extents[f_time_horizon][f_number_arcs]);
 for( Index t = 0; t < f_time_horizon; ++t ) {
 for( Index g = 0; g < f_number_arcs; ++g ) {
 v_secondary_spinning_reserve[ t ][ g ].set_type( ColVariable::kNonNegative );

  }
 }
  add_static_variable ( v_secondary_spinning_reserve, "sr" );
 // AR |= HasVar;
} // end( HydroUnitMaintenance::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void HydroUnitMaintenance::generate_abstract_constraints( Configuration *stcc ) {

 // maximum power output according to primary-secondary reserves constraints

 ///////////////////////////////////////////////////////////////////////////
 //                 Add the constraint of maintenance
 ///////////////////////////////////////////////////////////////////////////
    
  if( MaxPowerPrimarySecondary_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( MaxPowerPrimarySecondary_Const.empty());

   MaxPowerPrimarySecondary_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][f_number_arcs] );
  }
 for( Index arc = 0; arc < f_number_arcs; ++arc ) {

  for( Index t = 0; t < f_time_horizon; ++t ) {
    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_active_power[t][arc], 1.0 );
    if ( !v_primary_rho.empty() ) {
     linear_function->add_variable( &v_primary_spinning_reserve[t][arc], 1.0 );
    }
    if ( !v_secondary_rho.empty() ) {
     linear_function->add_variable( &v_secondary_spinning_reserve[t][arc], 1.0 );
    }
    
    // Here we add the v_maintenance (= 1 - \xi) variable
    if (!v_maximum_power.empty() ) {
        linear_function->add_variable( &v_maintenance[t][arc], - v_maximum_power[t][arc] );
    } else {
        linear_function->add_variable( &v_maintenance[t][arc], - v_maximum_flow[t][arc] );
    }
      
//    if (!v_minimum_power.empty() ) {
//     MaxPowerPrimarySecondary_Const[t][arc].set_lhs( v_minimum_power[t][arc] );
//    } else {
//     MaxPowerPrimarySecondary_Const[t][arc].set_lhs( 0.0 );
//
//    }

    MaxPowerPrimarySecondary_Const[t][arc].set_rhs( 0.0 );

    MaxPowerPrimarySecondary_Const[t][arc].set_function( linear_function );
   }
  }
  add_static_constraint( MaxPowerPrimarySecondary_Const, "MaxPowerPrimarySecondary" );
    
    
    
//////////////////////////////////////////////////////////////////////////////////////
//                     Add v_maintenance = z (= 0 for eg.) constraint               //
//////////////////////////////////////////////////////////////////////////////////////
    
    
    SplittingVar_Const.resize
               ( boost::multi_array< FRowConstraint, 2 >::
                 extent_gen()[f_time_horizon][f_number_arcs] );


    for( Index arc = 0; arc < f_number_arcs; ++arc ) {

      for( Index t = 0; t < f_time_horizon; ++t ) {

        auto linear_function = new LinearFunction();

        linear_function->add_variable( &v_maintenance[t][arc], 1.0 );

        SplittingVar_Const[t][arc].set_lhs( 0.0 );

        SplittingVar_Const[t][arc].set_rhs( 0.0 );

      }
    }
    
    add_static_constraint( SplittingVar_Const, "SplittingVar" );
  
//////////////////////////////////////////////////////////////////////////////////////
        
        

 // minimum power output according to primary-secondary reserves constraints

 if( MinPowerPrimarySecondary_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( MinPowerPrimarySecondary_Const.empty());

  MinPowerPrimarySecondary_Const.resize
          ( boost::multi_array< FRowConstraint, 2 >::
            extent_gen()[f_time_horizon][f_number_arcs] );
 }
 for( Index arc = 0; arc < f_number_arcs; ++arc ) {

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t][arc], 1.0 );
   if ( !v_primary_rho.empty() ) {
    linear_function->add_variable( &v_primary_spinning_reserve[t][arc], -1.0 );
   }
   if ( !v_secondary_rho.empty() ) {
    linear_function->add_variable( &v_secondary_spinning_reserve[t][arc], -1.0 );
   }
   if (!v_minimum_power.empty() ) {
    MinPowerPrimarySecondary_Const[t][arc].set_lhs( v_minimum_power[t][arc] );
   } else {
    MinPowerPrimarySecondary_Const[t][arc].set_lhs( 0.0 );

   }
   if (!v_maximum_power.empty() ) {
    MinPowerPrimarySecondary_Const[t][arc].set_rhs( v_maximum_power[t][arc] );
   } else {
    MinPowerPrimarySecondary_Const[t][arc].set_rhs( v_linear_term[0] * v_maximum_flow[t][arc] );

   }
   MinPowerPrimarySecondary_Const[t][arc].set_function( linear_function );
  }
 }

 add_static_constraint( MinPowerPrimarySecondary_Const, "MinPowerPrimarySecondary");

 // power output relation with to primary reserves constraints
 if( v_maximum_flow[0][0] > 0 ) {

  if( !v_primary_rho.empty()) {
   if( ActivePowerPrimary_Const.size() != f_time_horizon ) {
    // this should only happen once
    assert( ActivePowerPrimary_Const.empty());

    ActivePowerPrimary_Const.resize
            ( boost::multi_array< FRowConstraint, 2 >::
              extent_gen()[f_time_horizon][f_number_arcs] );
   }
   for( Index arc = 0; arc < f_number_arcs; ++arc ) {

    for( Index t = 0; t < f_time_horizon; ++t ) {

     auto linear_function = new LinearFunction();

     linear_function->add_variable( &v_active_power[t][arc], v_primary_rho[t][arc] );

     linear_function->add_variable( &v_primary_spinning_reserve[t][arc], -1.0 );
     ActivePowerPrimary_Const[t][arc].set_lhs( 0.0 );
     ActivePowerPrimary_Const[t][arc].set_rhs( Inf< double >());
     ActivePowerPrimary_Const[t][arc].set_function( linear_function );
    }

   }
   add_static_constraint( ActivePowerPrimary_Const, "ActivePowerPrimary" );

  }
 }
 // power output relation with to secondary reserves constraints
 if( v_maximum_flow[0][0] > 0 ) {
  if( !v_secondary_rho.empty()) {
   if( ActivePowerSecondary_Const.size() != f_time_horizon ) {
    // this should only happen once
    assert( ActivePowerSecondary_Const.empty());

    ActivePowerSecondary_Const.resize
            ( boost::multi_array< FRowConstraint, 2 >::
              extent_gen()[f_time_horizon][f_number_arcs] );
   }
   for( Index arc = 0; arc < f_number_arcs; ++arc ) {
    for( Index t = 0; t < f_time_horizon; ++t ) {
     auto linear_function = new LinearFunction();
     linear_function->add_variable( &v_active_power[t][arc], v_secondary_rho[t][arc] );
     linear_function->add_variable( &v_secondary_spinning_reserve[t][arc], -1.0 );
     ActivePowerSecondary_Const[t][arc].set_lhs( 0.0 );
     ActivePowerSecondary_Const[t][arc].set_rhs( Inf< double >());
     ActivePowerSecondary_Const[t][arc].set_function( linear_function );
    }
   }
   add_static_constraint( ActivePowerSecondary_Const, "ActivePowerSecondary" );
  }
 }
 // primary reserves constraints for pumps
 if( v_maximum_flow[0][0] <= 0 ) {
  if( PrimaryPumps_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( PrimaryPumps_Const.empty());

   PrimaryPumps_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][f_number_arcs] );
  }
  for( Index arc = 0; arc < f_number_arcs; ++arc ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {


     auto linear_function = new LinearFunction();

     linear_function->add_variable( &v_primary_spinning_reserve[t][arc], 1.0 );
     PrimaryPumps_Const[t][arc].set_both( 0.0 );
     PrimaryPumps_Const[t][arc].set_function( linear_function );
    }

  }

  add_static_constraint( PrimaryPumps_Const, "PrimaryPumps");
 }
 // secondary reserves constraints for pumps
 if( v_maximum_flow[0][0] <= 0 ) {

  if( SecondaryPumps_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( SecondaryPumps_Const.empty());

   SecondaryPumps_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][f_number_arcs] );
  }
  for( Index arc = 0; arc < f_number_arcs; ++arc ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {

     auto linear_function = new LinearFunction();

     linear_function->add_variable( &v_secondary_spinning_reserve[t][arc], 1.0 );
     SecondaryPumps_Const[t][arc].set_both( 0.0 );
     SecondaryPumps_Const[t][arc].set_function( linear_function );
    }

  }

  add_static_constraint( SecondaryPumps_Const, "SecondaryPumps");
 }

 // flow to active power function constraints for pumps
 if( v_maximum_flow[0][0] <= 0 ) {

  if( FlowActivePowerPumps_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( FlowActivePowerPumps_Const.empty());

   FlowActivePowerPumps_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][f_number_arcs] );
  }
  for( Index arc = 0; arc < f_number_arcs; ++arc ) {

  for( Index t = 0; t < f_time_horizon; ++t ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_active_power[t][arc], 1.0 );
    linear_function->add_variable( &v_flow_rate[t][arc], -v_linear_term[arc] );
    FlowActivePowerPumps_Const[t][arc].set_both( 0.0 );
    FlowActivePowerPumps_Const[t][arc].set_function( linear_function );

   }
  }
  add_static_constraint( FlowActivePowerPumps_Const, "FlowActivePowerPumps");
 }

 // flow to active power function constraints for turbines
 if (v_maximum_flow[0][0] > 0 ) {
  if (f_number_arcs > 0 ) {
  if( FlowActivePowerTurbines_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( FlowActivePowerTurbines_Const.empty());

   if( f_total_number_pieces == f_number_arcs ) {

    FlowActivePowerTurbines_Const.resize
            ( boost::multi_array< FRowConstraint, 3 >::
              extent_gen()[f_time_horizon][f_number_arcs][1] );
   } else if( f_total_number_pieces > f_number_arcs ) {
    FlowActivePowerTurbines_Const.resize
            ( boost::multi_array< FRowConstraint, 3 >::
              extent_gen()[f_time_horizon][1][f_total_number_pieces] );
   }
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {

   Index piece = 0;
   Index end = 0;
   for( Index arc = 0; arc < f_number_arcs; ++arc ) {
    if( !v_number_pieces.empty()) {
     end += v_number_pieces[arc];
    }
    if( f_total_number_pieces == f_number_arcs ) {

     auto linear_function = new LinearFunction();

     linear_function->add_variable( &v_active_power[t][arc], 1.0 );
     if( !v_linear_term.empty()) {
      linear_function
              ->add_variable( &v_flow_rate[t][arc], -v_linear_term[arc] );
     } else {
      linear_function->add_variable( &v_flow_rate[t][arc], 0.0 );

     }
     if( !v_const_term.empty()) {
      FlowActivePowerTurbines_Const[t][arc][0]
              .set_rhs( v_const_term[arc] );
     } else {
      FlowActivePowerTurbines_Const[t][arc][0].set_rhs( 0.0 );

     }
     FlowActivePowerTurbines_Const[t][arc][0].set_lhs( -Inf< double >());
     FlowActivePowerTurbines_Const[t][arc][0]
             .set_function( linear_function );

    } else if( f_total_number_pieces > f_number_arcs ) {

     for( ; piece < end; ++piece ) {
      auto linear_function = new LinearFunction();

      linear_function->add_variable( &v_active_power[t][arc], 1.0 );
      if( !v_linear_term.empty()) {
       linear_function
               ->add_variable( &v_flow_rate[t][arc], -v_linear_term[piece] );
      } else {
       linear_function->add_variable( &v_flow_rate[t][arc], 0.0 );
      }
      if( !v_const_term.empty()) {
       FlowActivePowerTurbines_Const[t][0][piece]
               .set_rhs( v_const_term[piece] );
      } else {
       FlowActivePowerTurbines_Const[t][0][piece].set_rhs( 0.0 );
      }
      FlowActivePowerTurbines_Const[t][0][piece]
              .set_lhs( -Inf< double >());
      FlowActivePowerTurbines_Const[t][0][piece]
              .set_function( linear_function );
     }
    }
   }
  }
  add_static_constraint( FlowActivePowerTurbines_Const, "FlowActivePowerTurbines" );
 }
 }
 // flow rate bounds constraints
 {
  if( FlowRateBounds_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( FlowRateBounds_Const.empty());

   FlowRateBounds_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][f_number_arcs] );
  }
  for( Index arc = 0; arc < f_number_arcs; ++arc ) {

   for( Index t = 0; t < f_time_horizon; ++t ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_flow_rate[t][arc], 1.0 );
    FlowRateBounds_Const[t][arc].set_lhs( 0.0 );
    FlowRateBounds_Const[t][arc].set_rhs( v_maximum_flow[t][arc] );
    FlowRateBounds_Const[t][arc].set_function( linear_function );

   }
  }

  add_static_constraint( FlowRateBounds_Const, "FlowRateBounds" );
 }


 // ram-up constraints
 if( !v_delta_ramp_up.empty() ) {

  if( RampUp_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( RampUp_Const.empty());

   RampUp_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][f_number_arcs] );
  }
  // Initial condition

  for( Index arc = 0; arc < f_number_arcs; ++arc ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_flow_rate[0][arc], 1.0 );

   RampUp_Const[0][arc].set_lhs( -Inf< double >());
   RampUp_Const[0][arc].set_rhs( v_delta_ramp_up[0][arc] + v_initial_flow_rate[arc] );
   RampUp_Const[0][arc].set_function( linear_function );

  }

  for( Index t = 1; t < f_time_horizon; ++t ) {
   for( Index arc = 0; arc < f_number_arcs; ++arc ) {
    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_flow_rate[t][arc], 1.0 );
    linear_function->add_variable( &v_flow_rate[t - 1 ][arc], -1.0 );

    RampUp_Const[t][arc].set_lhs( -Inf< double >());
    RampUp_Const[t][arc].set_rhs( v_delta_ramp_up[t][arc] );
    RampUp_Const[t][arc].set_function( linear_function );

   }
  }

  add_static_constraint( RampUp_Const,"RampUp");
 }


 // ram-down constraints
 if( !v_delta_ramp_down.empty()) {

  if( RampDown_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( RampDown_Const.empty());

   RampDown_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][f_number_arcs] );
  }

  // Initial condition

  for( Index arc = 0; arc < f_number_arcs; ++arc ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_flow_rate[0][arc], 1.0 );

   RampDown_Const[0][arc].set_lhs( v_initial_flow_rate[arc] - v_delta_ramp_down[0][arc] );
   RampDown_Const[0][arc].set_rhs( Inf< double >());
   RampDown_Const[0][arc].set_function( linear_function );

  }

  for( Index arc = 0; arc < f_number_arcs; ++arc ) {
  for( Index t = 1; t < f_time_horizon; ++t ) {
    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_flow_rate[t-1][arc], 1.0 );
    linear_function->add_variable( &v_flow_rate[t][arc], -1.0 );

    RampDown_Const[t][arc].set_lhs( -Inf< double >());
    RampDown_Const[t][arc].set_rhs( v_delta_ramp_down[t][arc] );
    RampDown_Const[t][arc].set_function( linear_function );

   }
  }

  add_static_constraint( RampDown_Const, "RampDown");

 }

 // final volumes fo each reservoir constraints

 if( FinalVolumeReservoir_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( FinalVolumeReservoir_Const.empty());
  FinalVolumeReservoir_Const.resize
          ( boost::multi_array< FRowConstraint, 2 >::
            extent_gen()[f_time_horizon][ f_number_reservoirs ] );
 }
 Index end = 0;

 for( Index n = 0; n < f_number_reservoirs; ++n ) {

  auto l_f = new LinearFunction();

  l_f->add_variable( &v_volumetric[n][0], 1.0 );

  end += n;

  for( Index l = 0; l <= end; ++l ) {

   if( !v_start_arc.empty() && !v_end_arc.empty()) {

     if( v_start_arc[l] == n ) {

      if( !v_uphill_delay.empty()) {

       l_f->add_variable( &v_flow_rate[0 - v_uphill_delay[l]][l], 1.0 );

       } else {

       l_f->add_variable( &v_flow_rate[0][l], 1.0 );

         }
   } else  {

      if( !v_downhill_delay.empty()) {
       if ( v_downhill_delay[l] == 0 ){
       l_f->add_variable( &v_flow_rate[0][l], -1.0 );
       }

      } else {

       l_f->add_variable( &v_flow_rate[0][l], -1.0 );

      }
   }

   } else {

     l_f->add_variable( &v_flow_rate[0][l], 1.0 );

    }

  }
  FinalVolumeReservoir_Const[0][n].set_both( v_initial_volumetric[n] + v_inflows[n][0] );
  FinalVolumeReservoir_Const[0][n].set_function( l_f );

  for( Index t = 1, constraint_index = 1; t < f_time_horizon;
       ++t, ++constraint_index  ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_volumetric[n][t], 1.0 );
   linear_function->add_variable( &v_volumetric[n][t-1], -1.0 );

   for( Index l = 0; l <= end; ++l ) {

    if( !v_start_arc.empty() && !v_end_arc.empty()) {


      if( v_start_arc[l] == n ) {

      if( !v_uphill_delay.empty()) {

       linear_function->add_variable( &v_flow_rate[t - v_uphill_delay[l]][l], 1.0 );

      } else {

       linear_function->add_variable( &v_flow_rate[t][l], 1.0 );

      }
     } else  {

      if( !v_downhill_delay.empty()  ) {

        linear_function->add_variable( &v_flow_rate[t - v_downhill_delay[l]][l], -1.0 );

      } else {

       linear_function->add_variable( &v_flow_rate[t][l], -1.0 );
      }
      }

    } else {

     linear_function->add_variable( &v_flow_rate[t][l], 1.0 );

    }
   }
   FinalVolumeReservoir_Const[constraint_index][n].set_both(  v_inflows[n][constraint_index] );
   FinalVolumeReservoir_Const[constraint_index][n].set_function( linear_function );
  }
  }

 add_static_constraint( FinalVolumeReservoir_Const, "FinalVolumeReservoir");



 // volumetric bounds constraints

 if( VolumetricBounds_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( VolumetricBounds_Const.empty());

  VolumetricBounds_Const.resize
          ( boost::multi_array< FRowConstraint, 2 >::
            extent_gen()[f_number_reservoirs][f_time_horizon] );
 }
 for( Index node = 0; node < f_number_reservoirs; ++node ) {
  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();
   linear_function->add_variable( &v_volumetric[node][t], 1.0 );

   VolumetricBounds_Const[node][t].set_lhs( v_minimum_volumetric[node][t] );
   VolumetricBounds_Const[node][t].set_rhs( v_maximum_volumetric[node][t] );
   VolumetricBounds_Const[node][t].set_function( linear_function );

  }
 }

 add_static_constraint( VolumetricBounds_Const, "VolumetricBounds");
 //AR |= HasCst;
} // end( HydroUnitMaintenance::generate_abstract_constraints )


/////////////////////////////////////////////////////////////////////////////
//                Generate objective for HydroUnitMaintenance              //
/////////////////////////////////////////////////////////////////////////////


void HydroUnitMaintenance::generate_objective( Configuration *objc ) {

    // Initialize objective function

    auto linear_function = new LinearFunction();

    for( Index arc = 0; arc < f_number_arcs; ++arc ) {

      for( Index t = 0; t < f_time_horizon; ++t ) {

            auto cost = 1.0;

            linear_function->add_variable( & v_active_power[t][arc], cost );

        }
    }

    objective.set_function( linear_function );
    objective.set_sense( Objective::eMin );
    objective.set_Block( this );

}


/*--------------------------------------------------------------------------*/
/*-------- METHODS FOR LOADING, PRINTING & SAVING THE HydroUnitMaintenance -------*/
/*--------------------------------------------------------------------------*/

void HydroUnitMaintenance::serialize( netCDF::NcGroup & group ) const {

 UnitBlock::serialize( group );

 auto dim_time_horizon = group.addDim( "TimeHorizon", f_time_horizon );
 auto NumberIntervals = group.getDim( "NumberIntervals" );


 auto dim_total_number_pieces = group.addDim( "TotalNumberPieces",
                                              f_total_number_pieces );

 auto dim_number_reservoirs = group.addDim( "NumberReservoirs",
                                            f_number_reservoirs
                                            ? f_number_reservoirs : 1 );

 auto dim_number_arcs = group.addDim( "NumberArcs",
                                      f_number_arcs ? f_number_arcs : 1 );


 ::serialize( group, "NumberPieces", netCDF::NcUint(),
              dim_number_arcs, v_number_pieces, true );

 ::serialize( group, "StartLine", netCDF::NcUint(),
              dim_number_reservoirs, v_start_arc, false );

 ::serialize( group, "EndLine", netCDF::NcUint(),
              dim_number_reservoirs, v_end_arc, false );


 if( !v_minimum_flow.empty() ) {

  ::serialize( group, "MinFlow", netCDF::NcDouble(),
               { NumberIntervals, dim_number_arcs },
               v_minimum_flow, true );
 }

 if( !v_maximum_flow.empty() ) {

  ::serialize( group, "MaxFlow", netCDF::NcDouble(),
               { NumberIntervals, dim_number_arcs },
               v_maximum_flow, true );
 }

 if( !v_minimum_volumetric.empty() ) {

  ::serialize( group, "MinVolumetric", netCDF::NcDouble(),
               { dim_number_reservoirs, NumberIntervals },
               v_minimum_volumetric, true );
 }

 if( !v_maximum_volumetric.empty() ) {

  ::serialize( group, "MaxVolumetric", netCDF::NcDouble(),
               { dim_number_reservoirs, NumberIntervals },
               v_maximum_volumetric, true );
 }

 ::serialize( group, "Inflows", netCDF::NcDouble(),
              { dim_number_reservoirs, dim_time_horizon },
              v_inflows, false );

 ::serialize( group, "MinPower", netCDF::NcDouble(),
              { NumberIntervals, dim_number_arcs },
              v_minimum_power, true );

 ::serialize( group, "MaxPower", netCDF::NcDouble(),
              { NumberIntervals, dim_number_arcs },
              v_maximum_power, true );

 ::serialize( group, "DeltaRampUp", netCDF::NcDouble(),
              { NumberIntervals, dim_number_arcs },
              v_delta_ramp_up, true );

 ::serialize( group, "DeltaRampDown", netCDF::NcDouble(),
              { NumberIntervals, dim_number_arcs },
              v_delta_ramp_down, true );

 if( !v_primary_rho.empty() ) {
  ::serialize( group, "PrimaryRho", netCDF::NcDouble(),
               { NumberIntervals, dim_number_arcs },
               v_primary_rho, true );
 }

 if( !v_secondary_rho.empty() ) {

  ::serialize( group, "SecondaryRho", netCDF::NcDouble(),
               { NumberIntervals, dim_number_arcs },
               v_secondary_rho, true );
 }

 ::serialize( group, "LinearTerm", netCDF::NcDouble(),
              dim_total_number_pieces, v_linear_term, false );

 ::serialize( group, "ConstantTerm", netCDF::NcDouble(),
              dim_total_number_pieces, v_const_term, false );


 ::serialize( group, "InertiaPower", netCDF::NcDouble(),
              { NumberIntervals, dim_number_arcs },
              v_inertia_power, true );

 ::serialize( group, "InitialFlowRate", netCDF::NcDouble(),
              dim_number_arcs, v_initial_flow_rate, false );

 ::serialize( group, "InitialVolumetric", netCDF::NcDouble(),
              dim_number_reservoirs, v_initial_volumetric, false );

 ::serialize( group, "UphillFlow", netCDF::NcUint(),
              dim_number_arcs, v_uphill_delay, true );

 ::serialize( group, "DownhillFlow", netCDF::NcUint(),
              dim_number_arcs, v_downhill_delay, true );
}  // end( HydroUnitMaintenance::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/
void
HydroUnitMaintenance::set_inflow( std::vector< double >::const_iterator values,
                            Block::Subset && subset,
                            const bool ordered,
                            c_ModParam issuePMod,
                            c_ModParam issueAMod ) {
 if( subset.empty() ) {
  return;
 }

 if( v_inflows.empty() ) {
  if( std::all_of( values,
                   values + subset.size(),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_inflows.resize( boost::extents[ f_number_reservoirs ][ f_time_horizon ] );
 }

 // If nothing changes, return
 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_inflows.size() ) {
   throw ( std::invalid_argument( "invalid value in subset" ) );
  }
  if( *( v_inflows.data() + i ) != *( values++ ) ) {
   identical = false;
  }
 }
 if( identical ) {
  return;
 }
 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  for( auto i : subset ) {
   Index t = i % f_time_horizon;
   Index r = i / f_time_horizon;
   v_inflows[ r ][ t ] = *( values++ );
   // *( v_inflows.data() + i ) = *( values++ );
  }

  /*if( AR & HasCst ) {
   // Change the abstract representation

   for( auto i : subset ) {
    Index t = i % f_time_horizon;
    Index r = i / f_time_horizon;

    if( t == 0 ) {
     FinalVolumeReservoir_Const[ t ][ r ]
      .set_both( v_initial_volumetric[ r ] + v_inflows[ r ][ t ], issueAMod );
    } else {
     FinalVolumeReservoir_Const[ t ][ r ]
      .set_both( v_inflows[ r ][ t ], issueAMod );
    }
   }
  }
 }
*/
 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered ) {
   std::sort( subset.begin(), subset.end() );
  }

  Block::add_Modification(
   std::make_shared< HydroUnitMaintenanceSbstMod >( this,
                                              HydroUnitMaintenanceMod::eSetInf,
                                              std::move( subset ) ),
   Observer::par2chnl( issuePMod ) );
 }
}

void
HydroUnitMaintenance::set_inflow( std::vector< double >::const_iterator values,
                            Block::Range rng,
                            c_ModParam issuePMod,
                            c_ModParam issueAMod ) {
 rng.second = std::min( rng.second,
                        get_time_horizon() * get_number_reservoirs() );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_inflows.empty() ) {
  if( std::all_of( values,
                   values + ( rng.second - rng.first ),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_inflows.resize( boost::extents[ f_number_reservoirs ][ f_time_horizon ] );
 }

 // If nothing changes, return
 if( std::equal( values,
                 values + ( rng.second - rng.first ),
                 v_inflows.data() + rng.first ) ) {
  return;
 }
 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values,
             values + ( rng.second - rng.first ),
             v_inflows.data() + rng.first );

  /*if( AR & HasCst ) {
   // Change the abstract representation

   for( Index i = rng.first; i < rng.second; ++i ) {
    Index t = i % f_time_horizon;
    Index r = i / f_time_horizon;

    if( t == 0 ) {
     FinalVolumeReservoir_Const[ t ][ r ]
      .set_both( v_initial_volumetric[ r ] + v_inflows[ r ][ t ], issueAMod );
    } else {
     FinalVolumeReservoir_Const[ t ][ r ]
      .set_both( v_inflows[ r ][ t ], issueAMod );
    }
   }
  }
 }*/

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification(
   std::make_shared< HydroUnitMaintenanceRngdMod >( this,
                                              HydroUnitMaintenanceMod::eSetInf,
                                              rng ),
   Observer::par2chnl( issuePMod ) );
 }
}

void
HydroUnitMaintenance::set_inertia_power( std::vector< double >::const_iterator values,
                                   Subset && subset,
                                   const bool ordered,
                                   c_ModParam issuePMod,
                                   c_ModParam issueAMod ) {
 if( subset.empty() ) {
  return;
 }

 if( v_inertia_power.empty() ) {
  if( std::all_of( values,
                   values + subset.size(),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_inertia_power.resize( boost::extents[ f_time_horizon ][ f_number_arcs ] );
 }

 // If nothing changes, return
 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_inertia_power.size() ) {
   throw ( std::invalid_argument( "invalid value in subset" ) );
  }
  if( *( v_inertia_power.data() + i ) != *( values++ ) ) {
   identical = false;
  }
 }
 if( identical ) {
  return;
 }
 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  for( auto i : subset ) {
   Index a = i % f_number_arcs;
   Index t = i / f_number_arcs;
   v_inertia_power[ t ][ a ] = *( values++ );
  }

  /*if( AR & HasCst ) {
   // Change the abstract representation
   // FIXME: v_inertia_power is not used
  }
 }*/

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered ) {
   std::sort( subset.begin(), subset.end() );
  }

  Block::add_Modification(
   std::make_shared< HydroUnitMaintenanceSbstMod >( this,
                                              HydroUnitMaintenanceMod::eSetInerP,
                                              std::move( subset ) ),
   Observer::par2chnl( issuePMod ) );
 }
}

void
HydroUnitMaintenance::set_inertia_power( std::vector< double >::const_iterator values,
                                   Block::Range rng,
                                   c_ModParam issuePMod,
                                   c_ModParam issueAMod ) {
 rng.second = std::min( rng.second, f_number_arcs * get_time_horizon() );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_inertia_power.empty() ) {
  if( std::all_of( values,
                   values + ( rng.second - rng.first ),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_inertia_power.resize( boost::extents[ f_time_horizon ][ f_number_arcs ] );
 }

 // If nothing changes, return
 if( std::equal( values,
                 values + ( rng.second - rng.first ),
                 v_inertia_power.data() + rng.first ) ) {
  return;
 }
 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values,
             values + ( rng.second - rng.first ),
             v_inertia_power.data() + rng.first );

  /*if( AR & HasCst ) {
   // Change the abstract representation
   // FIXME: v_inertia_power is not used

  }
 }*/

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification(
   std::make_shared< HydroUnitMaintenanceRngdMod >( this,
                                              HydroUnitMaintenanceMod::eSetInerP,
                                              rng ),
   Observer::par2chnl( issuePMod ) );
 }
}

void
HydroUnitMaintenance::set_initial_volumetric(
 std::vector< double >::const_iterator values,
 Block::Subset && subset,
 const bool ordered,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 if( subset.empty() ) {
  return;
 }

 if( v_initial_volumetric.empty() ) {
  if( std::all_of( values,
                   values + subset.size(),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  Index max_index = *max_element( std::begin( subset ), std::end( subset ) );
  v_initial_volumetric.assign( max_index, 0 );
 }

 // If nothing changes, return
 bool identical = true;
 auto temp_values = values;
 for( auto i : subset ) {
  if( i >= v_initial_volumetric.size() ) {
   throw ( std::invalid_argument( "invalid value in subset" ) );
  }
  if( v_initial_volumetric[ i ] != *( temp_values++ ) ) {
   identical = false;
  }
 }
 if( identical ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  temp_values = values;
  for( auto i : subset ) {
   v_initial_volumetric[ i ] = *( temp_values++ );
  }

  /*if( not_dry_run( issueAMod ) && AR & HasObj ) {
   // Change the abstract representation
   for( auto i : subset ) {
    Index t = i % f_time_horizon;
    Index r = i / f_time_horizon;

    if( t == 0 ) {
     FinalVolumeReservoir_Const[ t ][ r ]
      .set_both( v_initial_volumetric[ r ] + v_inflows[ r ][ t ], issueAMod );
    }
   }
  }
 }
 */
     
 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered ) {
   std::sort( subset.begin(), subset.end() );
  }
  Block::add_Modification(
   std::make_shared< HydroUnitMaintenanceSbstMod >( this,
                                              HydroUnitMaintenanceMod::eSetInitV,
                                              std::move( subset ) ),
   Observer::par2chnl( issuePMod ) );
 }
}

void
HydroUnitMaintenance::set_initial_volumetric(
 std::vector< double >::const_iterator values,
 Block::Range rng,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 rng.second = std::min( rng.second, f_time_horizon );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_initial_volumetric.empty() ) {
  if( std::all_of( values,
                   values + ( rng.second - rng.first ),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  Index max_index = rng.second;
  v_initial_volumetric.assign( max_index, 0 );
 }

 // If nothing changes, return
 if( std::equal( values,
                 values + ( rng.second - rng.first ),
                 v_initial_volumetric.begin() + rng.first ) ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values,
             values + ( rng.second - rng.first ),
             v_initial_volumetric.begin() + rng.first );

  /*if( not_dry_run( issueAMod ) && AR & HasCst ) {
   // Change the abstract representation
   for( Index i = rng.first; i < rng.second; ++i ) {
    Index t = i % f_time_horizon;
    Index r = i / f_time_horizon;

    if( t == 0 ) {
     FinalVolumeReservoir_Const[ t ][ r ]
      .set_both( v_initial_volumetric[ r ] + v_inflows[ r ][ t ], issueAMod );
    }
   }
  }
 }*/

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification(
   std::make_shared< HydroUnitMaintenanceRngdMod >( this,
                                              HydroUnitMaintenanceMod::eSetInitV,
                                              rng ),
   Observer::par2chnl( issuePMod ) );
 }
}

template< typename T >
void HydroUnitMaintenance::transpose( boost::multi_array< T, 2 > & a ) {
 long rows = a.shape()[ 0 ];
 long cols = a.shape()[ 1 ];
 if( rows > 1 && cols == 1 ) {
  // The vector must be transposed
  boost::array< typename boost::multi_array< T, 2 >::index, 2 > dims = { { 1, rows } };
  a.reshape( dims );
 }
}

void HydroUnitMaintenance::decompress_array( boost::multi_array< double, 2 > & a ) {

 if (a.empty()) {
  return;
 }
 boost::multi_array< double, 2 > temp = a;
 a.resize( boost::extents[ f_time_horizon ][ f_number_arcs ] );

 if( a.shape()[ 1 ] == 1 ) {
  for( Index t = 0; t < f_time_horizon; ++t ) {
   for( Index g = 0; g < f_number_arcs; ++g ) {
    a[ t ][ g ] = temp[ 0 ][ g ];
   }
  }

 } else if( a.shape()[ 0 ] < f_time_horizon ) {
  for( Index g = 0; g < f_number_arcs; ++g ) {
   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[ i ];
    }
    for( ; j < sup; ++j ) {
     a[ j ][ g ] = temp[ i ][ g ];
    }
   }
  }
 }
}

void HydroUnitMaintenance::decompress_vol( boost::multi_array< double, 2 > & a ) {

 if (a.empty()) {
  return;
 }
 boost::multi_array< double, 2 > temp = a;
 a.resize( boost::extents[ f_number_reservoirs ][ f_time_horizon ] );

 if( a.shape()[ 0 ] == 1 ) {

  for( Index n = 0; n < f_number_reservoirs; ++n ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {
    a[ n ][ t ] = temp[ n ][ 0 ];
   }
  }

 } else if( a.shape()[ 1 ] < f_time_horizon ) {
  for( Index n = 0; n < f_number_reservoirs; ++n ) {
   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[ i ];
    }
    for( ; j < sup; ++j ) {
     a[ n ][ j ] = temp[ n ][ i ];
    }
   }
  }
 }
}
/*--------------------------------------------------------------------------*/
/*------------------- End File HydroUnitMaintenance.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
