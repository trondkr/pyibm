import datetime as datetime

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2010, 6, 15)
__modified__ = datetime.datetime(2010, 6, 15)
__version__  = "1.0"
__status__   = "Development"


def showInfo():
    print """
    Introduction
    ------------
    This is an Individual-Based Model developed for Atlantic cod (Gadus morhua).
    The model was written by Trond Kristiansen.
    First version of this model was written in Fortran code as part of Trond Kristiansen
    PhD thesis (see: http://www.assembla.com/code/ibm/subversion/nodes). Starting in
    the summer of 2009, the model wasre-written and strongly improved as part of
    the ARCWARM project where we wanted to use the IBM to model larval cod growth and survival
    for a variety of time periods and geographic locations. The model was therefore improved by
    using netCDF input files that specify timeunits as "Days since (1948.1.1)". So, days
    is the de-facto standard timeunit for this model and all input files should use that
    format to make sure all interpolations in time and space are valid.

    THE IBM
    --------
    This IBM simulate early life history of larval cod based on environmental
    factors and underlying processes that we believe are important. Each process
    is formulated and parameterized from laboratory experiments conducted on larval
    cod (Fiksen and MacKenzie 2002). The combination of these processes describes
    mechanistically how larval cod encounters, capture, and ingest
    food and use the energy for growth. The main purpose is to evaluate model
    performance in a well-studied semi-natural situation.

    Sub-models govern the interaction between the environmental forcing
    (temperature, turbulence, and zooplankton density) and the dynamic state
    variables weight, stomach content, and length. These states are updated
    once every hour (Fig. 2). A key element is the mechanistic foraging
    sub-model combined with a stomach as a state variable. Ingested mass or
    energy results from the sequential processes; encounter, approach, and the
    probability of successful capture (Fiksen & MacKenzie 2002). The foraging
    processes are iterated for all prey types and size categories (Fig. 2).
    Holling's disc equation (Holling 1966) ensures that the total time spent
    handling prey reduces available search time. Stomach fullness and body mass
    determine whether growth is only temperature-dependent, or also food limited.

    The model is referenced and described in detail in these refs:
    --------------------------------------------------------------

    Kristiansen T., Fiksen O. & Folkvord A. (2007).
    Modelling feeding, growth, and habitat selection in larval Atlantic cod
    (Gadus morhua): observations and model predictions in a macrocosm
    environment. Can. J. Fish. Aquat. Sci., 64, 136-151.

    Kristiansen T., Jorgensen C., Lough R.G., Vikeboe F. & Fiksen O. (2009).
    Modeling rule-based behavior: habitat selection and the growth-survival
    trade-off in larval cod. Behav. Ecol., 20, 490-500.

    Kristiansen T., Lough R.G., Werner F.E., Broughton E.A. &
    Buckley L.J. (2009). Individual-based modeling of feeding
    ecology and prey selection of larval cod on Georges Bank.
    Mar. Ecol. Progr. Ser., 376, 227-243.

    Kristiansen T. (2007). Modeling early life history of cod.
    Department of Biology, Ph.D., 248.

    """
def showZooplanktonInfo():
    print """\nUSE ZOOPLANKTON MODEL:
    The option 'useZooplanktonModel' is turned on (init function in ibm.py). This means
    that we will calculate the zooplankton density for each time step and depth layer
    using the small and large phytoplankton production values together with
    temperature from the ESM model. The actual calculation is done in
    the Fortran routine calculateZooplankton found in the bioenergetics module.

    This approach assumes that the order of the variables
    (in routine getStationData in module IOnetcdf.py) is exactly:
    varlist=['temp','salt','nh4sm','nh4lg','no3sm','no3lg','chla','taux','tauy']
    If this is not correct, change the order in init function and re-run.\n
    """
