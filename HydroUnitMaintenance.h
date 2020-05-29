/*--------------------------------------------------------------------------*/
/*------------------------- File HydroUnitMaintenance.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class HydroUnitMaintenance, which derives from UnitBlock
 * [see UnitBlock.h], in order to define a "reasonably standard" hydro unit
 * of a Unit Commitment Problem.
 *
 * \version 0.11
 *
 * \date 19 - 05 - 2020
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __HydroUnitMaintenance
#define __HydroUnitMaintenance
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ColVariable.h"
#include "FRowConstraint.h"
#include "DQuadFunction.h"
#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// Namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it {


class HydroUnitMaintenance : public UnitBlock {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:
/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC TYPES OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/// constructor, takes the father and the time horizon
/** Constructor of HydroUnitMaintenance, taking possibly a pointer of its
 * father Block.
 */

 explicit HydroUnitMaintenance( Block * f_block = nullptr , Index t = 0):
         UnitBlock( f_block , t ) {}

/*--------------------------------------------------------------------------*/

/// destructor of HydroUnitMaintenance

 ~HydroUnitMaintenance() override = default;



 void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// generate the abstract variables of the HydroUnitMaintenance
/** The HydroUnitMaintenance class has five boost::multi_array< ColVariable, 2 >
 *  variables where the first four are:
 *
 *  - the primary spinning reserve variables;
 *
 *  - the secondary spinning reserve variables;
 *
 *  - the active power variables;
 *
 *  - the flow rate variables
 *
 *   Each of the boost::multi_array< ColVariable, 2 > has as first dimension
 *   the time horizon and as second dimension the number of arcs(or generators
 *   which as returned get_number_generators()). The last
 *   boost::multi_array< ColVariable, 2 > variable is:
 *
 *  - the volumetric variables
 *
 *  Where its' first dimension is the number of reservoirs and the second
 *  dimension is the time horizon.
 *
 *  All of these variables are optional,and it is also possible to restrict
 *  which of the subsets are generated with the parameter stvv. If stvv is not
 *  nullptr and it is a SimpleConfiguration<int>, or if
 *  f_BlockConfig->f_static_variables_Configuration is not nullptr and it is a
 *  SimpleConfiguration<int>, then the f_value (an int) indicates whether each
 *  of the optional variables should be created. If the Configuration is not
 *  available, the default value is taken to be 0.
 * */        
    
/* Here we need to add the maintenance variables who are the same dimension as active power variables */
    
    
 void generate_abstract_variables( Configuration *stvv ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the static constraint of the HydroUnit
/** This method generates the static constraint of the HydroUnitMaintenance.
 * In order to describe a hydro generating unit system, it will be convenient
 * to see a cascading system as a graph. Let \f$ \mathcal{N}^{hy}\f$ be the
 * set of reservoirs (nodes) and \f$ \mathcal{L}^{hy}\f$ be the set of arcs
 * connecting these reservoirs respectively. Attached to each
 * \f$ l \in \mathcal{L}^{hy}\f$ are one or several plants(turbines or pumps).
 * This system is described on a discrete time horizon as dictated by the
 * UnitBlock interface. In this description we indicate it with
 * \f$ \mathcal{T}=\{ 0, \dots , \mathcal{|T|} - 1\} \f$. Each reservoir
 * \f$ n \in \mathcal{N}^{hy}\f$ has a continuous volumetric variables
 * \f$ v^{hy}_{n,t}\f$ in \f$ m^3 \f$ for \f$ t \in \mathcal{T}\f$ with
 * associated lower and upper bounds \f$ V^{hy,mn}_{n,t}\f$,
 * \f$ V^{hy,mx}_{n,t}\f$ and inflows \f$ A_{n,t}\f$ in \f$ m^3 /s \f$. The
 * uphill and downhill flow rate are defined as \f$ \tau^{up} \f$ and
 * \f$ \tau^{dn}\f$ respectively. For each arc \f$ l \in \mathcal{L}^{hy}\f$
 * in each time \f$ t \in \mathcal{T}\f$ the continuous flow rate variable
 * \f$ f_{l,t} \f$ in \f$ m^3 /s \f$ and ramping conditions
 * \f$ \Delta^{up}_{l,t} \f$ and \f$ \Delta^{dn}_{l,t} \f$ in
 * \f$ (m^3 /s)/h \f$ are disposed. The flow rate variable will be subject to
 * bounds \f$ F^{mn}_{l,t} \f$ and \f$ F^{mx}_{l,t} \f$ and it's assumed
 * moreover given a cutting plane model describing power as a function of flow
 * rate as below:
 *   \f[
 *      p^{ac}_{t,l}(f) := min \{ P_l + \rho^{hy}_{l}f_{t,l}\}
 *   \f]
 * where \f$ P_l\f$ and \f$ \rho^{hy}_{l} \f$ are considered as the constant
 * and linear multipliers of the linear function(flow-to-active-power)
 * respectively. Power generated by the hydro unit in each time and for each
 * arc \f$ p^{ac}_{t,l}, p^{pr}_{t,l}, p^{sc}_{t,l} \f$ in MW will be subject
 * to bounds \f$ P^{mn}_{t,l} \f$ and \f$ P^{mx}_{t,l} \f$ respectively.
 * Besides, we emphasize that reserve requirements are specified in order to
 * be symmetrically available to increase or decrease power injected into the
 * grid. For some of the constraints we will need to distinguish between pumps
 * and turbines. The distinction is made by considering the set of feasible
 * flow rates. Whenever \f$ [ F^{mn}_{l,t} , F^{mx}_{l,t}] \subseteq R_- \f$
 * for each arc and each time, the unit is considered a pump, and whenever
 * \f$ [ F^{mn}_{l,t} , F^{mx}_{l,t}] \subseteq R_+ \f$ the unit is considered
 * a turbine. Any possible mixed situation can be accounted for by
 * artificially splitting the unit into “two units”, which should be done at
 * the data processing stage (see deserialize() comments). With above
 * description the mathematical constraint of hydro unit may present as below:
 *
 * - maximum and minimum power output constraints according to primary and
 *   secondary spinning reserves are are presented in (1)-(2). Each of them
 *   is a boost::multi_array<FRowConstraint, 2>; with two dimensions which are
 *   f_time_horizon, and f_number_arcs entries, where the entry
 *   t = 0, ...,f_time_horizon - 1 and the entry l = 0, ...,f_number_arcs - 1
 *   being the maximum and minimum power output value according to the primary
 *   and the secondary spinning reserves at time t and arc l. these ensure the
 *   maximum(or minimum) amount of energy that unit can produce(or use) when
 *   it is on(or off).
 *
 *   \f[
 *
 *      p^{ac}_{t,l} + p^{pr}_{t,l} + p^{sc}_{t,l} \leq P^{mx}_{t,l}
 *          \quad t \in \mathcal{T}, l \in \mathcal{L}^{hy}          \quad (1)
 *
 *   \f]
 *
 *   \f[
 *
 *     P^{mn}_{t,l} \leq p^{ac}_{t,l} - p^{pr}_{t,l} - p^{sc}_{t,l}
 *         \quad t \in \mathcal{T}, l \in \mathcal{L}^{hy}           \quad (2)
 *
 *   \f]
 *
 * - primary and secondary spinning reserves relation with active power at
 *   each time and for each turbine: the same as inequalities(1)-(2), the
 *   inequalities(3)-(4) ensure that maximum amount of primary and secondary
 *   spinning reserve in the problem. Each of them is a
 *   boost::multi_array<FRowConstraint, 2>; with two dimensions which are
 *   f_time_horizon, and f_number_arcs entries, where
 *   t = 0, ...,f_time_horizon - 1 and l = 0, ...,f_number_arcs - 1
 *
 *   \f[
 *
 *      p^{pr}_{t,l} \leq \rho^{pr}_{t,l}p^{ac}_{t,l} \quad t \in \mathcal{T},
 *        l \in \mathcal{L}^{hy} \quad with
 *        \quad  [ F^{mn}_{l,t} , F^{mx}_{l,t}] \subseteq R_+       \quad (3)
 *
 *   \f]
 *
 *   \f[
 *
 *      p^{sc}_{t,l} \leq \rho^{sc}_{t,l}p^{ac}_{t,l} \quad t \in \mathcal{T},
 *        l \in \mathcal{L}^{hy} \quad with
 *        \quad  [ F^{mn}_{l,t} , F^{mx}_{l,t}] \subseteq R_+       \quad (4)
 *
 *   \f]
 *  where \f$ \rho^{pr}_{t,l} \f$ and \f$ \rho^{sc}_{t,l}\f$ are the maximum
 *  possible fraction of active power at each time and each arc that can be
 *  used as primary and secondary reserve respectively.
 *
 * - primary and secondary spinning reserves at each time and for each pump:
 *   these equalities(5)-(6) ensure that the primary and secondary spinning
 *   reserve value for each pump is equal to zero. Each of them is a
 *   boost::multi_array<FRowConstraint, 2>; with two dimensions which are
 *   f_time_horizon, and f_number_arcs entries, where
 *   t = 0, ...,f_time_horizon - 1 and l = 0, ...,f_number_arcs - 1
 *
 *   \f[
 *
 *      p^{pr}_{t,l} = 0 \quad t \in \mathcal{T},
 *        l \in \mathcal{L}^{hy} \quad with
 *        \quad  [ F^{mn}_{l,t} , F^{mx}_{l,t}] \subseteq R_-      \quad (5)
 *
 *   \f]
 *
 *   \f[
 *
 *      p^{sc}_{t,l} = 0 \quad t \in \mathcal{T},
 *        l \in \mathcal{L}^{hy} \quad with
 *        \quad  [ F^{mn}_{l,t} , F^{mx}_{l,t}] \subseteq R_-     \quad (6)
 *
 *   \f]
 *
 * - flow to active power function at each time and for each pump: this
 *   equality(7) gives the active power relation with flow rate for each
 *   pump at time t. This is a boost::multi_array<FRowConstraint, 2>; with two
 *   dimensions which are f_time_horizon, and f_number_arcs entries, where
 *   t = 0, ...,f_time_horizon - 1 and l = 0, ...,f_number_arcs - 1
 *
 *   \f[
 *
 *      p^{ac}_{t,l} = \rho^{hy}_{l}f_{t,l} \quad t \in \mathcal{T},
 *        l \in \mathcal{L}^{hy} \quad with
 *        \quad  [ F^{mn}_{l,t} , F^{mx}_{l,t}] \subseteq R_-       \quad (7)
 *
 *   \f]
 *
 * - flow-to-active-power function at each time and for each turbine ;
 *
 *   \f[
 *
 *      p^{ac}_{t,l} \leq P_j + \rho^{hy}_{j}f_{t,l} \quad j \in \mathcal{J}_l
 *        \quad t \in \mathcal{T}, l \in \mathcal{L}^{hy} \quad with
 *        \quad  [ F^{mn}_{l,t} , F^{mx}_{l,t}] \subseteq R_+       \quad (8)
 *
 *   \f]
 *
 * - flow rate variable bounds: This inequality(9) indicates upper and lower
 *   bound of flow rate at time t and for ach arc l, thus that is a
 *   boost::multi_array<FRowConstraint, 2>; with two dimensions which are
 *   f_time_horizon, and f_number_arcs entries, where
 *   t = 0, ...,f_time_horizon - 1 and l = 0, ...,f_number_arcs - 1
 *
 *   \f[
 *
 *      f_{t,l} \in [ F^{mn}_{t,l} , F^{mx}_{t,l}]  \quad t \in \mathcal{T},
 *           l \in \mathcal{L}^{hy}                             \quad (9)
 *
 *   \f]
 *
 * - ramp-up and ramp-down constraints: These inequality(10)-(11) indicate
 *   ramp-up and ramp-down constraints at time t and for ach arc l, so each of
 *   them is a boost::multi_array<FRowConstraint, 2>; with two dimensions
 *   which are f_time_horizon, and f_number_arcs entries, where
 *   t = 0, ...,f_time_horizon - 1 and l = 0, ...,f_number_arcs - 1
 *
 *   \f[
 *
 *      f_{t,l} - f_{t-1,l} \leq \Delta^{up}_{t,l} \quad t \in \mathcal{T},
 *           l \in \mathcal{L}^{hy}                             \quad (10)
 *
 *   \f]
 *   \f[
 *
 *      f_{t-1,l} - f_{t,l} \leq \Delta^{dn}_{t,l} \quad t \in \mathcal{T},
 *           l \in \mathcal{L}^{hy}                             \quad (11)
 *
 *   \f]
 *
 * - final volumes of each reservoir constraints: this equality(12) gives the
 *   final volumes of each reservoir r at time t. This is a
 *   boost::multi_array<FRowConstraint, 2>; with two dimensions which are
 *   f_number_reservoirs, and f_time_horizon entries, where
 *   n = 0, ...,f_number_reservoirs - 1 and t = 0, ...,f_time_horizon - 1
 *   \f[
 *
 *      v^{hy}_{n,t} = v^{hy}_{n,t-1} + A_{n,t} + (\sum_{l=(n',n) \in
 *      \mathcal{L}^{hy} } f_{t - \tau^{dn}_l , l } - \sum_{l=(n,n') \in
 *      \mathcal{L}^{hy} } f_{t - \tau^{up}_l , l }) \quad t \in \mathcal{T},
 *      \quad n \in \mathcal{N}^{hy}    \quad (12)
 *
 *   \f]
 *   where in each arc \f$ l=(n,n') \in \mathcal{L}^{hy} \f$, \f$ n \f$ and
 *   \f$ n' \f$ are supposed to be the start and the end point of that
 *   respectively.
 *
 * - final volumes variable bounds: This inequality(13) indicates upper and
 *   lower bound of volumetric variables of each reservoir for each time t,
 *   thus that is a boost::multi_array<FRowConstraint, 2>; with two dimensions
 *   which are f_number_reservoirs, and f_time_horizon entries, where
 *   n = 0, ...,f_number_reservoirs - 1 and t = 0, ...,f_time_horizon - 1
 *   \f[
 *
 *      v^{hy}_{n,t} \in [ V^{hy,mn}_{n,t} , V^{hy,mx}_{n,t}]  \quad
 *        n \in \mathcal{N}^{hy}, t \in \mathcal{T}            \quad (13)
 *
 *   \f]
 *
 */
 void generate_abstract_constraints( Configuration *stcc ) override;
    
    
 void generate_objective( Configuration *objc )  override;
/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE HydroUnitMaintenance -------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the HydroUnitMaintenance
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of hydro units
 * @{ */

/// returns the number of reservoirs
 Index get_number_reservoirs() const {
  return f_number_reservoirs? f_number_reservoirs : 1 ; }

/*--------------------------------------------------------------------------*/
 /// returns the number of arcs/generators
 Index get_number_generators() const override {
  return f_number_arcs ? f_number_arcs : 1;
 }

/*--------------------------------------------------------------------------*/
/// returns the vector of start arcs
/** Method for returning the vector of starting point of each arc. This vector
 * may have size of 1 (single hydro unit with just one arc between two
 * reservoirs) or the size of number of reservoirs, then there are two
 * possible cases:
 *
 *  - if f_number_reservoirs == 2, this vector has size of 1 which means there
 *    is one arc at the system (single hydro case).
 *
 *  - if f_number_reservoirs > 2, this vector have size of f_number_reservoirs
 *    and each element of the vectors gives starting point of each arc in the
 *    network.
 */
 const std::vector< Index > & get_start_arc() const { return v_start_arc; }

/*--------------------------------------------------------------------------*/
/// returns the vector of end arcs
/** Method for returning the vector of ending point of each arc. This
 *  vector may have size of 1 (single hydro unit with just one arc between two
 *  reservoirs) or the size of number of reservoirs, then there are two
 *  possible cases:
 *
 *  - if f_number_reservoirs == 2, this vector has size of 1 which means there
 *    is one arc at the system (single hydro case).
 *
 *  - if f_number_reservoirs > 2, this vector have size of f_number_reservoirs
 *    and each element of the vectors gives ending point of each arc in the
 *    network.
 */
 const std::vector< Index > & get_end_arc() const { return( v_end_arc ); }
 /*--------------------------------------------------------------------------*/

/// returns the matrix of inertia power
/** The returned value U = get_inertia_power() contains the contribution to
 *  inertia (basically, the constants to be multiplied by the active power
 *  variables returned by get_active_power()) of eac arcs (generators) at each
 *  time instants. There are three possible cases
 *
 * - if the matrix is empty, then the inertia power is 0;
 *
 * - if the matrix only has one row (i.e., the first dimension has size 1),
 *   then the inertia power for arc (generator) l is U[ 0 , l ] for all t
 *   which means that the second dimension has size get_number_arcs();
 *
 * - otherwise, the matrix has size the time horizon per
 *   get_number_arcs(), and U[ t , l ] contains the contribution to
 *   inertia power of arc(generator) l at time instant t. */
 double * get_inertia_power( Index generator )  override {
  return ( v_inertia_power.data() + generator * f_time_horizon );
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of minimum volumetric
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ n , t ] gives the minimum volumetric of the reservoir n at the time
 * instant t. This two-dimensional boost::multi_array<> M considers three
 * possible cases:
 *
 *  - if the boost::multi_array<> M is empty() then no minimum volumetric are
 *    defined, and there are no minimum volumetric constraints;
 *
 *  - if the boost::multi_array<> M has only one column which in this case the
 *    boost::multi_array<> M is a (transpose of) vector with size
 *    get_number_reservoirs(). Each element of M[ n , 0 ] gives the minimum
 *    volumetric of the each reservoir n for all time instant t;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_number_reservoirs() row where each row must have size of
 *    get_time_horizon() and each element of M[ n , t ] gives the minimum
 *    volumetric of reservoir n at time instant t. */
 const boost::multi_array< double , 2 > & get_minimum_volumetric() const {
  return ( v_minimum_volumetric );
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of maximum volumetric
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ n , t ] gives the maximum volumetric of the reservoir n at the time
 * instant t. This two-dimensional boost::multi_array<> M considers three
 * possible cases:
 *
 *  - if the boost::multi_array<> M is empty() then no maximum volumetric are
 *    defined, and there are no minimum volumetric constraints;
 *
 *  - if the boost::multi_array<> M has only one column which in this case the
 *    boost::multi_array<> M is a (transpose of) vector with size
 *    get_number_reservoirs(). Each element of M[ n , 0 ] gives the maximum
 *    volumetric of the each reservoir n for all time instant t;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_number_reservoirs() row where each row must have size of
 *    get_time_horizon() and each element of M[ n , t ] gives the maximum
 *    volumetric of reservoir n at time instant t. */
 const boost::multi_array< double , 2 > & get_maximum_volumetric() const {
  return ( v_maximum_volumetric );
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of inflows
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ n , t ] gives the inflows of the reservoir n at the time
 * instant t. This two-dimensional boost::multi_array<> M considers two
 * possible cases:
 *
 *  - if the boost::multi_array<> M is empty() then no inflows are defined;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_number_reservoirs() row where each row must have size of
 *    get_time_horizon() and each element of M[ n , t ] represents the inflows
 *    of reservoir n at time instant t. */
 const boost::multi_array< double , 2 > & get_inflows() const {
  return ( v_inflows );
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of minimum power
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ t , i ] gives the minimum power at each time t associated with unit(arc)
 * i. This two-dimensional boost::multi_array<> M considers three possible
 * cases:
 *
 *  - if the boost::multi_array<> M is empty() then no minimum power are
 *    defined, and there are no minimum power constraints;
 *
 *  - if the boost::multi_array<> M has only one row which in this case the
 *    boost::multi_array<> M is a vector with size get_number_arcs(). Each
 *    element of M[ 0 , i ] gives the minimum power for all time instant t of
 *    each unit i;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_time_horizon() rows and get_number_arcs() columns and each element
 *    of M[ t , i ] gives the minimum power at time t and unit i. */
 const boost::multi_array< double , 2 > & get_minimum_power() const {
  return ( v_minimum_power );
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of maximum power
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ t , i ] gives the maximum power at each time t associated with unit(arc)
 * i. This two-dimensional boost::multi_array<> M considers three possible
 * cases:
 *
 *  - if the boost::multi_array<> M is empty() then no maximum power are
 *    defined, and there are no maximum power constraints;
 *
 *  - if the boost::multi_array<> M has only one row which in this case the
 *    boost::multi_array<> M is a vector with size get_number_arcs(). Each
 *    element of M[ 0 , i ] gives the maximum power for all time instant t of
 *    each unit i;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_time_horizon() rows and get_number_arcs() columns and each element
 *    of M[ t , i ] gives the maximum power at time t and unit i. */
 const boost::multi_array< double , 2 > & get_maximum_power() const {
  return ( v_maximum_power );
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of minimum flow
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ t , i ] gives the minimum flow at each time t associated with unit(arc)
 * i. This two-dimensional boost::multi_array<> M considers three possible
 * cases:
 *
 *  - if the boost::multi_array<> M is empty() then no minimum flow are
 *    defined, and there are no minimum flow constraints;
 *
 *  - if the boost::multi_array<> M has only one row which in this case the
 *    boost::multi_array<> M is a vector with size get_number_arcs(). Each
 *    element of M[ 0 , i ] gives the minimum flow for all time instant t of
 *    each unit i;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_time_horizon() rows and get_number_arcs() columns and each element
 *    of M[ t , i ] gives the minimum flow at time t and unit i. */
 const boost::multi_array< double , 2 > & get_minimum_flow() const {
  return ( v_minimum_flow );
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of maximum flow
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ t , i ] gives the maximum flow at each time t associated with unit(arc)
 * i. This two-dimensional boost::multi_array<> M considers three possible
 * cases:
 *
 *  - if the boost::multi_array<> M is empty() then no maximum flow are
 *    defined, and there are no maximum flow constraints;
 *
 *  - if the boost::multi_array<> M has only one row which in this case the
 *    boost::multi_array<> M is a vector with size get_number_arcs(). Each
 *    element of M[ 0 , i ] gives the maximum flow for all time instant t of
 *    each unit i;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_time_horizon() rows and get_number_arcs() columns and each element
 *    of M[ t , i ] gives the maximum flow at time t and unit i. */
 const boost::multi_array< double , 2 > & get_maximum_flow() const {
  return ( v_maximum_flow );
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of delta ramp up
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ t , i ] gives the delta ramp up at each time t associated with unit(arc)
 * i. This two-dimensional boost::multi_array<> M considers three possible
 * cases:
 *
 *  - if the boost::multi_array<> M is empty() then no ramping constraints;
 *
 *  - if the boost::multi_array<> M has only one row which in this case the
 *    boost::multi_array<> M is a vector with size get_number_arcs(). Each
 *    element of M[ 0 , i ] gives the delta ramp up value for all time instant
 *    t of each unit i;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_time_horizon() rows and get_number_arcs() columns and each element
 *    of M[ t , i ] gives the delta ramp up value at time t and unit i. */
 const boost::multi_array< double , 2 > & get_delta_ramp_up() const {
  return ( v_delta_ramp_up);
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of delta ramp down
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ t , i ] gives the delta ramp down at each time t associated with unit
 * (arc) i. This two-dimensional boost::multi_array<> M considers three
 * possible cases:
 *
 *  - if the boost::multi_array<> M is empty() then no ramping constraints;
 *
 *  - if the boost::multi_array<> M has only one row which in this case the
 *    boost::multi_array<> M is a vector with size get_number_arcs(). Each
 *    element of M[ 0 , i ] gives the delta ramp down value for all time
 *    instant t of each unit i;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_time_horizon() rows and get_number_arcs() columns and each element
 *    of M[ t , i ] gives the delta ramp down value at time t and unit i. */
 const boost::multi_array< double , 2 > & get_delta_ramp_down() const {
  return ( v_delta_ramp_down);
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of primary rho
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ t , i ] gives the primary rho at each time t associated with unit
 * (arc) i. This two-dimensional boost::multi_array<> M considers three
 * possible cases:
 *
 *  - if the boost::multi_array<> M is empty() then no primary reserve
 *    constraints;
 *
 *  - if the boost::multi_array<> M has only one row which in this case the
 *    boost::multi_array<> M is a vector with size get_number_arcs(). Each
 *    element of M[ 0 , i ] gives the primary rho value for all time instant t
 *    of each unit i;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_time_horizon() rows and get_number_arcs() columns and each element
 *    of M[ t , i ] gives the primary rho value at time t and unit i. */
 const boost::multi_array< double , 2 > & get_primary_rho() const {
  return ( v_primary_rho);
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of secondary rho
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ t , i ] gives the secondary rho at each time t associated with unit
 * (arc) i. This two-dimensional boost::multi_array<> M considers three
 * possible cases:
 *
 *  - if the boost::multi_array<> M is empty() then no secondary reserve
 *    constraints;
 *
 *  - if the boost::multi_array<> M has only one row which in this case the
 *    boost::multi_array<> M is a vector with size get_number_arcs(). Each
 *    element of M[ 0 , i ] gives the secondary rho value for all time instant
 *    t of each unit i;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_time_horizon() rows and get_number_arcs() columns and each element
 *    of M[ t , i ] gives the secondary rho value at time t and unit i. */
 const boost::multi_array< double , 2 > & get_secondary_rho() const {
  return ( v_secondary_rho);
 }
 /*--------------------------------------------------------------------------*/
 /// returns the vector of number pieces
/** The returned vector contains the number of pieces for each unit (arc) i.
 * There are three possible cases:
 *
 * - if the vector is empty, then the number of pieces of each unit is 0;
 *
 * - if the vector has only one element, then V[ 0 ] presents the number
 *   of pieces for all units;
 *
 * - otherwise, the vector must have size get_number_arcs() and each element
 *   of V[ i ] represents the number of pieces of each unit i. */
 const std::vector< Index > & get_number_pieces() const {
  return ( v_number_pieces);
 }
 /*--------------------------------------------------------------------------*/
 /// returns the vector of constant term
/** The returned vector contains the constant term value of each available
 * pieces h in the set of {0, ..., TotalNumberPieces}(see NumberPieces
 * deserialize() comments). There are three possible cases:
 *
 * - if the vector is empty, then the constant term value is 0;
 *
 * - if the vector has only one element, then V[ 0 ] presents the const term
 *   value for all the pieces;
 *
 * - otherwise, the returned V is a std::vector < double > and
 *   V.sized == TotalNumberPieces and each element of V[ h ] represents the
 *   const term value of each piece h. */
 const std::vector< double > & get_const_term() const {
  return ( v_const_term);
 }
/*--------------------------------------------------------------------------*/
 /// returns the vector of linear term
/** The returned vector contains the linear term value of each available
 * pieces h in the set of {0, ..., TotalNumberPieces}(see NumberPieces
 * deserialize() comments). There are three possible cases:
 *
 * - if the vector is empty, then the linear term value is 0;
 *
 * - if the vector has only one element, then V[ 0 ] presents the linear term
 *   value for all the pieces;
 *
 * - otherwise, the returned V is a std::vector < double > and
 *   V.sized == TotalNumberPieces and each element of V[ h ] represents the
 *   linear term value of each piece h. */
 const std::vector< double > & get_linear_term() const {
  return ( v_linear_term);
 }
/*--------------------------------------------------------------------------*/
 /// returns the vector of uphill delay
/** The returned vector contains the uphill delay for each unit (arc) i.
 * There are three possible cases:
 *
 * - if the vector is empty, then the uphill delay for each unit is 0;
 *
 * - if the vector has only one element, then V[ 0 ] presents the uphill delay
 *   for all units;
 *
 * - otherwise, the vector must have size get_number_arcs() and each element
 *   of V[ i ] represents the uphill delay for each unit i. */
 const std::vector< Index > & get_uphill_delay() const {
  return ( v_uphill_delay);
 }
/*--------------------------------------------------------------------------*/
 /// returns the vector of downhill delay
/** The returned vector contains the downhill delay for each unit (arc) i.
 * There are three possible cases:
 *
 * - if the vector is empty, then the downhill delay for each unit is 0;
 *
 * - if the vector has only one element, then V[ 0 ] presents the downhill
 *   delay for all units;
 *
 * - otherwise, the vector must have size get_number_arcs() and each element
 *   of V[ i ] represents the downhill delay for each unit i. */
 const std::vector< Index > & get_downhill_delay() const {
  return ( v_downhill_delay);
 }
/*--------------------------------------------------------------------------*/
 /// returns the vector of initial volumetric
/** The returned vector contains the initial volumetric for each reservoir n.
 * There are three possible cases:
 *
 * - if the vector is empty, then the initial volumetric for each reservoir is
 *   0;
 *
 * - if the vector has only one element, then V[ 0 ] presents the initial
 *   volumetric for all reservoirs;
 *
 * - otherwise, the vector must have size get_number_reservoirs() and each
 *   element of V[ i ] represents the initial volumetric for each reservoir n.
 *   */
 const std::vector< double > & get_initial_volumetric() const {
  return ( v_initial_volumetric);
 }
/*--------------------------------------------------------------------------*/
 /// returns the vector of initial flow rate
/** The returned vector contains the initial flow rate for each unit (arc) i.
 * There are three possible cases:
 *
 * - if the vector is empty, then the initial flow rate for each unit is 0;
 *
 * - if the vector has only one element, then V[ 0 ] presents the initial flow
 *   rate for all units;
 *
 * - otherwise, the vector must have size get_number_arcs() and each element
 *   of V[ i ] represents the initial flow rate for each unit i. */
 const std::vector< double > & get_initial_flow_rate() const {
  return ( v_initial_flow_rate);
 }
/**@} ----------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE Variable OF THE HydroUnitMaintenancefget --------*/
/*--------------------------------------------------------------------------*/

/** @name Reading the Variable of the HydroUnitMaintenance
 *
 * These methods allow to read the five groups of Variable that any
 * HydroUnitMaintenance in principle has (although some may not):
 *
 *  - the volumetric variables;
 *
 *  - the flow rate variables;
 *
 *  - the active power variables;
 *
 *  - the primary spinning reserve variables;
 *
 *  - the secondary spinning reserve variables;
 *
 * All these five groups of variables are (if not empty)
 * boost::multi_array< ColVariable , 2 > where the volumetric variables with
 * the first dimension number reservoirs and the second dimension time horizon
 * and all the rest with first dimension time horizon and second dimension
 * number of arcs (generators).
 * @{ */

/// returns the matrix of volumetric variables
/** The returned boost::multi_array< ColVariable , 2 >, say V, contains the
 * volumetric variables and is indexed over the dimensions time horizon and
 * number of reservoirs. There are two possible cases:
 *
 * - if V is empty(), then these variables are not defined;
 *
 * - otherwise, V must have f_time_horizon rows and f_number_reservoir columns
 *   and M[ t , n ] is the volumetric variable for time step t of reservoir n.
 *   */

 ColVariable * get_volumetric( Index reservoir  ) {
  return ( v_volumetric.data() + reservoir * f_time_horizon );
 }

/*--------------------------------------------------------------------------*/

 /// returns the volume of the given reservoir at the given time
 /** Returns a pointer to the ColVariable representing the volume of the given
  * \p reservoir at the given \p time.
  *
  * @param reservoir The index of the reservoir whose volume is desired.
  *
  * @param time The time at which the volume is desired.
  *
  * @return A pointer to the ColVariable representing the volume of the given
  *         \p reservoir at the given \p time.
  */
 ColVariable * get_volume( Index reservoir , Index time ) {
  return ( v_volumetric.data() + reservoir * f_time_horizon + time );
 }

/*--------------------------------------------------------------------------*/
/// returns the matrix of flow rate variables
/** The returned boost::multi_array< ColVariable , 2 >, say F, contains the
 * flow rate variables and is indexed over the dimensions time horizon and
 * number of arcs (generators). There are two possible cases:
 *
 * - if F is empty(), then these variables are not defined;
 *
 * - otherwise, F must have f_time_horizon rows and f_number_arcs columns, and
 *   M[ t , a ] is the flow rate variable for time step t of arc (generator)
 *   a. */

 ColVariable * get_flow_rate( Index arc)  {
  return ( v_flow_rate.data() + arc * f_time_horizon);
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of active_power variables
 ColVariable * get_active_power( Index generator  ) override {
  return ( v_active_power.data() + generator * f_time_horizon );
 }
/*--------------------------------------------------------------------------*/
 /// returns the matrix of primary_spinning_reserve variables
 ColVariable * get_primary_spinning_reserve( Index generator )
 override {
  return ( v_primary_spinning_reserve.data() + generator * f_time_horizon );
 }
/*--------------------------------------------------------------------------*/
 /// returns the matrix of secondary_spinning_reserve variables
 ColVariable * get_secondary_spinning_reserve( Index generator )
 override {
  return ( v_secondary_spinning_reserve.data() + generator * f_time_horizon );
 }
/**@} ----------------------------------------------------------------------*/
/*------------------ METHODS FOR SAVING THE HydroUnitMaintenance------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the HydroUnitMaintenance
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * HydroUnitMaintenance. See HydroUnitMaintenance::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR INITIALIZING THE HydroUnitMaintenance --------------*/
/*--------------------------------------------------------------------------*/

/** @name Handling the data of the HydroUnitMaintenance
    @{ */

 void load( std::istream & input ) override {
  throw ( std::logic_error( "HydroUnitMaintenance::load() not implemented yet") );
 };

/**@} ----------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

 void set_inflow( std::vector< double >::const_iterator values,
                  Subset && subset,
                  bool ordered = false,
                  c_ModParam issuePMod = eNoBlck,
                  c_ModParam issueAMod = eNoBlck );

 void set_inflow( std::vector< double >::const_iterator values,
                  Range rng = Range( 0, Inf< Index >() ),
                  c_ModParam issuePMod = eNoBlck,
                  c_ModParam issueAMod = eNoBlck );

 void set_inertia_power( std::vector< double >::const_iterator values,
                         Subset && subset,
                         bool ordered = false,
                         c_ModParam issuePMod = eNoBlck,
                         c_ModParam issueAMod = eNoBlck );

 void set_inertia_power( std::vector< double >::const_iterator values,
                         Range rng = Range( 0, Inf< Index >() ),
                         c_ModParam issuePMod = eNoBlck,
                         c_ModParam issueAMod = eNoBlck );

 void set_initial_volumetric( std::vector< double >::const_iterator values,
                              Subset && subset,
                              bool ordered = false,
                              c_ModParam issuePMod = eNoBlck,
                              c_ModParam issueAMod = eNoBlck );

 void set_initial_volumetric( std::vector< double >::const_iterator values,
                              Range rng = Range( 0, Inf< Index >() ),
                              c_ModParam issuePMod = eNoBlck,
                              c_ModParam issueAMod = eNoBlck );



/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------data--------------------------------------*/
 /// The number of reservoirs(nodes) of the problem
 Index f_number_reservoirs;

 /// The number of connecting arcs which are connecting the reservoirs in
 /// cascading system
 Index f_number_arcs;

 /// The total number of pieces
 Index f_total_number_pieces;

 /// The vector of UphillDelay
 std::vector< Index > v_uphill_delay;

 /// The vector of DownhillDelay
 std::vector< Index > v_downhill_delay;

 /// The vector of starting arc
 std::vector< Index > v_start_arc;

 /// The vector of ending arcs
 std::vector< Index > v_end_arc;

 /// The vector of initial volumetric
 std::vector< double > v_initial_volumetric;

 /// The vector of initial flow rate
 std::vector< double > v_initial_flow_rate;

 /// The vector of NumberPieces
 std::vector< Index > v_number_pieces;

 /// The vector of LinearTerm
 std::vector< double > v_linear_term;

 /// The vector of ConstTerm
 std::vector< double > v_const_term;

 /// the matrix of inertia power of generators
 boost::multi_array< double , 2 > v_inertia_power;

 /// The matrix of MinVolumetric
 /** Indexed over the dimensions NumberReservoirs and NumberIntervals. */
 boost::multi_array< double, 2 > v_minimum_volumetric;

 /// The matrix of MaxVolumetric
 /** Indexed over the dimensions NumberReservoirs and NumberIntervals. */
 boost::multi_array< double, 2 > v_maximum_volumetric;

 /// The matrix of Inflows
 /** Indexed over the dimensions NumberReservoirs and NumberIntervals. */
 boost::multi_array< double, 2 > v_inflows;

 /// The matrix of MinPower
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double, 2 > v_minimum_power;

 /// The matrix of MaxPower
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double, 2 > v_maximum_power;

 /// The matrix of MinFlow
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double, 2 > v_minimum_flow;

 /// The matrix of MaxFlow
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double, 2 > v_maximum_flow;

 /// The matrix of DeltaRampUp
 /** Indexed over the dimensions NumberIntervals and NumberGenerators. */
 boost::multi_array< double, 2 > v_delta_ramp_up;

 /// The matrix of DeltaRampDown
 /** Indexed over the dimensions NumberIntervals and NumberGenerators. */
 boost::multi_array< double, 2 > v_delta_ramp_down;

 /// The matrix of PrimaryRho
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double, 2 > v_primary_rho;

 /// The matrix of SecondaryRho
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double, 2 > v_secondary_rho;

/*-----------------------------variables------------------------------------*/

 /// the matrix of volumetric variables
 boost::multi_array< ColVariable , 2 > v_volumetric;

 /// the matrix of flow rate variables
 boost::multi_array< ColVariable , 2 > v_flow_rate;

 /// the active power variables
 boost::multi_array< ColVariable , 2 > v_active_power;

 /// the primary spinning reserve variables
 boost::multi_array< ColVariable , 2 > v_primary_spinning_reserve;

 /// the secondary spinning reserve variables
 boost::multi_array< ColVariable , 2 > v_secondary_spinning_reserve;
/*----------------------------constraints-----------------------------------*/
 /// maximum power output according to primary-secondary reserves constraints
 boost::multi_array< FRowConstraint, 2 >  MaxPowerPrimarySecondary_Const;

 /// minimum power output according to primary-secondary reserves constraints
 boost::multi_array< FRowConstraint, 2 >  MinPowerPrimarySecondary_Const;

 /// power output relation with to primary reserves constraints
 boost::multi_array< FRowConstraint, 2 >  ActivePowerPrimary_Const;

 /// power output relation with to secondary reserves constraints
 boost::multi_array< FRowConstraint, 2 >  ActivePowerSecondary_Const;

 /// primary reserves constraints for pumps
 boost::multi_array< FRowConstraint, 2 >  PrimaryPumps_Const;

 /// secondary reserves constraints for pumps
 boost::multi_array< FRowConstraint, 2 >  SecondaryPumps_Const;

 /// flow to active power function constraints for pumps
 boost::multi_array< FRowConstraint, 2 >  FlowActivePowerPumps_Const;

 /// flow to active power function constraints for turbine
 boost::multi_array< FRowConstraint, 3 >  FlowActivePowerTurbines_Const;

 /// ramp-up constraints
 boost::multi_array< FRowConstraint, 2 >  RampUp_Const;

 /// ramp-down constraints
 boost::multi_array< FRowConstraint, 2 >  RampDown_Const;

 /// flow rate bounds constraints
 boost::multi_array< FRowConstraint, 2 >  FlowRateBounds_Const;

 /// final volumes fo each reservoir constraints
 boost::multi_array< FRowConstraint, 2 >  FinalVolumeReservoir_Const;

 /// volumetric bounds constraints
 boost::multi_array< FRowConstraint, 2 >  VolumetricBounds_Const;

 static void static_initialization() {
  register_method< HydroUnitMaintenance >( "HydroUnitMaintenance::set_inflow",
                                     &HydroUnitMaintenance::set_inflow,
                                     MS_dbl_sbst::args() );

  register_method< HydroUnitMaintenance >( "HydroUnitMaintenance::set_inflow",
                                     &HydroUnitMaintenance::set_inflow,
                                     MS_dbl_rngd::args() );

  register_method< HydroUnitMaintenance >( "HydroUnitMaintenance::set_inertia_power",
                                     &HydroUnitMaintenance::set_inertia_power,
                                     MS_dbl_sbst::args() );

  register_method< HydroUnitMaintenance >( "HydroUnitMaintenance::set_inertia_power",
                                     &HydroUnitMaintenance::set_inertia_power,
                                     MS_dbl_rngd::args() );

  register_method< HydroUnitMaintenance >( "HydroUnitMaintenance::set_initial_volumetric",
                                     &HydroUnitMaintenance::set_initial_volumetric,
                                     MS_dbl_sbst::args() );

  register_method< HydroUnitMaintenance >( "HydroUnitMaintenance::set_initial_volumetric",
                                     &HydroUnitMaintenance::set_initial_volumetric,
                                     MS_dbl_rngd::args() );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 /// FIXME this is also defined in UCBlock so let's put it in a common file
 /// Transposes a deserialized multiarray if needed.
 /**
  * Checks if the multiarray has one column and more than one rows. If so,
  * it transposes it.
  * This procedure is needed because some 2D matrices have the first dimension
  * optional, and the ::deserialize() method doesn't know that the only
  * dimension that is given is actually the second one.
  *
  * @tparam T The type of the boost::multi_array
  * @param a  A boost::multi_array that has been just deserialized
  */
 template< typename T >
 void transpose( boost::multi_array< T, 2 > & a );

 /// Decompress a multi_array using the change intervals
 void decompress_array( boost::multi_array< double, 2 > & a );

 /// Decompress a max/min volumetric multi_array using the change intervals
 void decompress_vol( boost::multi_array< double, 2 > & a );

/*--------------------------------------------------------------------------*/

};  // end( class( HydroUnitMaintenance ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS HydroUnitMaintenanceMod ------------------------*/
/*--------------------------------------------------------------------------*/

/// Derived class from Modification for modifications to a HydroUnitMaintenance
class HydroUnitMaintenanceMod : public Modification {

 public:

 /// Public enum for the types of HydroUnitMaintenanceMod
 enum HUB_mod_type {
  eSetInf = 0 ,    ///< Set inflow values
  eSetInerP    ,   ///< Set inertia power values
  eSetInitV        ///< Set initial volumetric values
 };

 /// Constructor, takes the HydroUnitMaintenance and the type
 HydroUnitMaintenanceMod( HydroUnitMaintenance * const fblock,
                      const int type )
  : f_Block( fblock ), f_type( type ) {}

 ///< Destructor, does nothing
 ~HydroUnitMaintenanceMod() override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block() const override { return ( f_Block ); }

 /// Accessor to the type of modification
 int type() { return ( f_type ); }

 protected:

 /// prints the HydroUnitMaintenanceMod
 void print( std::ostream & output ) const override {
  output << "HydroUnitMaintenanceMod[" << this << "]: ";
  switch( f_type ) {
   case ( eSetInf ):
    output << "set inflow values ";
    break;
   case ( eSetInerP ):
    output << "set inertia power values ";
    break;
   default:
    output << "set initial volumetric values ";
  }
 }

 HydroUnitMaintenance * f_Block{};
 ///< pointer to the Block to which the Modification refers

 int f_type; ///< type of modification
}; // end( class( HydroUnitMaintenanceMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS HydroUnitMaintenanceRngdMod ----------------------*/
/*--------------------------------------------------------------------------*/
/// derived from HydroUnitMaintenanceMod for "ranged" modifications
class HydroUnitMaintenanceRngdMod : public HydroUnitMaintenanceMod {

 public:

 /// constructor: takes the HydroUnitMaintenance, the type, and the range
 HydroUnitMaintenanceRngdMod( HydroUnitMaintenance * const fblock,
                          const int type,
                          Block::Range rng )
  : HydroUnitMaintenanceMod( fblock, type ), f_rng( rng ) {}

 /// destructor, does nothing
 ~HydroUnitMaintenanceRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng() { return( f_rng ); }

 protected:

 /// prints the HydroUnitMaintenanceRngdMod
 void print( std::ostream & output ) const override {
  HydroUnitMaintenanceMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
 }

 Block::Range f_rng; ///< the range
};  // end( class( HydroUnitMaintenanceRngdMod ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS HydroUnitMaintenanceSbstMod ---------------------*/
/*--------------------------------------------------------------------------*/

/// derived from HydroUnitMaintenanceMod for "subset" modifications
class HydroUnitMaintenanceSbstMod : public HydroUnitMaintenanceMod {

 public:

 /// constructor: takes the HydroUnitMaintenance, the type, and the subset
 HydroUnitMaintenanceSbstMod( HydroUnitMaintenance * const fblock,
                          const int type,
                          Block::Subset && nms )
  : HydroUnitMaintenanceMod( fblock, type ), f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 ~HydroUnitMaintenanceSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms() { return( f_nms ); }

 protected:

 /// prints the HydroUnitMaintenanceSbstMod
 void print( std::ostream &output ) const override {
  HydroUnitMaintenanceMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
 }

 Block::Subset f_nms; ///< the subset

};  // end( class( HydroUnitMaintenanceSbstMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* HydroUnitMaintenance.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File HydroUnitMaintenance.h -------------------------*/
/*--------------------------------------------------------------------------*/
