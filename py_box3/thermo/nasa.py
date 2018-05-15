from scipy.stats import variation
import numpy as np
import py_box3.constants as c
import warnings

class Nasa(object):
    """
    Contains the NASA polynomials and corresponding temperature ranges for a species.
    Information on NASA polynomials can be found here: http://combustion.berkeley.edu/gri_mech/data/nasa_plnm.html

    Attributes
    ----------
        T_low - float
            Lower temperature bound (K)
        T_mid - float
            Middle temperature where the transition between a_low and a_high occurs (K)
        T_high - float
            High temperature bound (K)
        a_low - (7,) ndarray 
            Holds the NASA coefficients used between T_low and T_mid.
        a_high - (7,) ndarray
            Holds the NASA coefficients used between T_mid and T_high.
        verbose - boolean
            Whether information should be printed as functions are called
    """
    def __init__(self, symbol = '', T_low = 0, T_mid= 0, T_high = 0, a_low = None, a_high = None, verbose = True):
        self.symbol = symbol
        self.T_low = T_low
        self.T_mid = T_mid
        self.T_high = T_high
        if a_low is None:
            a_low = np.array(7*[0.])
        self.a_low = a_low
        if a_high is None:
            a_high = np.array(7*[0.])
        self.a_high = a_high
        self.verbose = verbose


    def __str__(self):
        return str(self.__dict__)

    def __eq__(self, other):
        if all([all(self.a_high == other.a_high),
                all(self.a_low == other.a_low),
                self.T_high == other.T_high,
                self.T_mid == other.T_mid,
                self.T_low == other.T_low]):
            return True
        else:
            return False

    def _get_single_CpoR(self, T):
        """
        Calculates the heat capacity at constant pressure (i.e. Cp/R) for a single temperature.
        Parameters
        ----------
            T - float
                Temperature (K)
        Returns
        -------
            CpoR - float
                Heat capacity divided by molar gas constant (i.e. Cp/R)
        """
        T_arr = np.array([1., T, T ** 2, T ** 3, T ** 4, 0., 0.])
        if T < self.T_mid:
            if T < self.T_low:
                warnings.warn("Input temperature (%f) lower than T_low (%f)" % (T, self.T_low))
            return np.dot(T_arr, self.a_low)
        else:
            if T > self.T_high:
                warnings.warn("Warning. Input temperature (%f) higher than T_high (%f)" % (T, self.T_high))
            return np.dot(T_arr, self.a_high)


    def get_CpoR(self, T):
        """
        Calculates the heat capacity at constant pressure (i.e. Cp/R) given a temperature or a list of temperatures.
        Parameters
        ----------
            T - float or (N,) ndarray
                Temperature (K)
        Returns
        -------
            CpoR - float or (N,) ndarray
                Heat capacity divided by molar gas constant (i.e. Cp/R)
        """
        try:
            T_val = iter(T)
        except TypeError:
            #Single value T
            CpoR = self._get_single_CpoR(T)
        else:
            #List value T
            CpoR = np.array([0.]*len(T))
            for i, T_val in enumerate(T):
                CpoR[i] = self._get_single_CpoR(T_val)
        return CpoR

    def _get_single_HoRT(self, T, verbose = True):
        """
        Calculates the dimensionless enthalpy (i.e. H/RT) given a single temperature.
        Parameters
        ----------
            T - float
                Temperature (K)
        Returns
        -------
            HoRT - float
                Enthalpy divided by molar gas constant and temperature (i.e. H/RT)
        """
        T_arr = np.array([1., T/2., T ** 2 / 3., T ** 3 / 4., T ** 4 / 5., 1. / T, 0.])
        if T < self.T_mid:
            if T < self.T_low:
                warnings.warn("Input temperature (%f) lower than T_low (%f)" % (T, self.T_low))
            return np.dot(T_arr, self.a_low)
        else:
            if T > self.T_high:
                warnings.warn("Warning. Input temperature (%f) higher than T_high (%f)" % (T, self.T_high))
        return np.dot(T_arr, self.a_high)

    def get_HoRT(self, T):
        """
        Calculates the dimensionless enthalpy at constant pressure (i.e. H/RT) given a temperature or a list of temperatures.
        Parameters
        ----------
            T - float or (N,) ndarray
                Temperature (K)
        Returns
        -------
            HoRT - float or (N,) ndarray
                Enthalpy divided by molar gas constant and temperature (i.e. H/RT)
        """
        try:
            T_val = iter(T)
        except TypeError:
            #Single value T
            HoRT = self._get_single_HoRT(T)
        else:
            #List value T
            HoRT = np.array([0.]*len(T))
            for i, T_val in enumerate(T):
                HoRT[i] = self._get_single_HoRT(T_val)
        return HoRT


    def _get_single_SoR(self, T, verbose = True):
        """
        Calculates the dimensionless entropy (i.e. S/R) given a temperature.
        Parameters
        ----------
            T - float
                Temperature (K)
        Returns
        -------
            SoR - float
                Entropy divided by molar gas constant (i.e. S/R)
        """
        T_arr = np.array([np.log(T), T, T ** 2 / 2., T ** 3 / 3., T ** 4 / 4., 0., 1.])
        if T < self.T_mid:
            if T < self.T_low:
                warnings.warn("Input temperature (%f) lower than T_low (%f)" % (T, self.T_low))
            return np.dot(T_arr, self.a_low)
        else:
            if T > self.T_high:
                warnings.warn("Warning. Input temperature (%f) higher than T_high (%f)" % (T, self.T_high))
            return np.dot(T_arr, self.a_high)

    def get_SoR(self, T):
        """
        Calculates the dimensionless entropy at constant pressure (i.e. S/R) given a temperature or a list of temperatures.
        Parameters
        ----------
            T - float
                Temperature (K)
        Returns
        -------
            SoR - float
                Heat capacity divided by molar gas constant (i.e. S/R)
        """
        try:
            T_val = iter(T)
        except TypeError:
            #Single value T
            SoR = self._get_single_SoR(T)
        else:
            #List value T
            SoR = np.array([0.]*len(T))
            for i, T_val in enumerate(T):
                SoR[i] = self._get_single_SoR(T_val)
        return SoR

    def get_GoRT(self, T, verbose = True):
        """
        Calculates the dimensionless enthalpy (i.e. G/RT) given a temperature.
        Parameters
        ----------
            T - float or (N,) ndarray
                Temperature (K)
        Returns
        -------
            GoRT - float or (N,) ndarray
                Enthalpy divided by molar gas constant and temperature (i.e. G/RT)
        """
        HoRT = self.get_HoRT(T)
        SoR = self.get_SoR(T)
        return HoRT-SoR

    def plot_thermo(self, T_low = None, T_high = None, units = None):
        """
        Plots the heat capacity, enthalpy and entropy in the temperature range specified.
        Parameters
        ----------
            T_low - float
                Lower bound to plot temperature (K). If not specified then the T_low attribute for the species is used.
            T_high - float
                Higher bound to plot temperatures (K). If not specified then the T_high attribute for the species is used.
            units - string
                Controls the units used based on the molar gas constant (e.g. J/mol/K, kcal/mol/K, 'eV/K'). 
                If not specified then dimensionless units are used.
        """
        import matplotlib.pyplot as plt

        if T_low == None:
            T_low = self.T_low
        if T_high == None:
            T_high = self.T_high

        T = np.linspace(T_low, T_high)
        Cp = self.get_CpoR(T)
        H = self.get_HoRT(T)
        S = self.get_SoR(T)
        if units is not None:
            Cp = Cp * c.R(units)
            H = H * c.R(units) * T
            S = S * c.R(units)

        plt.figure()
        plt.subplot(311)
        plt.plot(T, Cp, 'r-')
        if units is None:
            plt.ylabel('Cp/R')
        else:
            plt.ylabel('Cp (%s)' % units)
        plt.xlabel('T (K)')
        plt.title('Plots for %s using NASA polynomials.' % self.symbol)

        plt.subplot(312)
        plt.plot(T, H, 'b-')
        if units is None:
            plt.ylabel('H/RT')
        else:
            plt.ylabel('H (%s)' % units.replace('/K', ''))
        plt.xlabel('T (K)')

        #Entropy graph
        plt.subplot(313)
        plt.plot(T, S, 'k-')
        if units is None:
            plt.ylabel('S/R')
        else:
            plt.ylabel('S (%s)' % units)
        plt.xlabel('T (K)')

    def fit_NASA(self, T, CpoR, HoRT0, SoR0, T0 = c.T0('K')):
        """
        Generate a NASA polynomial given dimensionless heat capacity as a function of temperature,
        dimensionless enthalpy of formation and dimensionless entropy of formation.
        Parameters
        ----------
            T - (N,) ndarray
                Temperatures (K) to fit the polynomial
            CpoR - (N,) ndarray
                Dimensionless heat capacities that correspond to T array
            HoRT0 - float
                Dimensionless enthalpy to be used as reference. Used to find a6. Value should correspond to T0
            SoR0 - float
                Dimensionless entropy to be used as reference. Used to find a7. Value should correspond to T0
            T0 - float
                Reference temperature used for fitting a6 and a7.
        """
        self.fit_CpoR(T = T, CpoR = CpoR)
        self._fit_HoRT(HoRT0 = HoRT0, T0 = T0)
        self._fit_SoR(SoR0 = SoR0, T0 = T0)

    def fit_CpoR(self, T, CpoR):
        """
        Fits parameters a1 - a6 using dimensionless heat capacity and temperature.
        Parameters
        ----------
            T - (N,) ndarray
                Temperatures (K) to fit the polynomial
            CpoR - (N,) ndarray
                Dimensionless heat capacities that correspond to T array
        """
        #If the Cp/R does not vary with temperature (occurs when no vibrational frequencies are listed)
        if (np.mean(CpoR) < 1e-6 and np.isnan(variation(CpoR))) or variation(CpoR) < 1e-3 or all(np.isnan(CpoR)):
           self.T_mid = T[int(len(T)/2)]
           self.a_low = np.array(7*[0.])
           self.a_high = np.array(7*[0.])
        else:
            max_R2 = -1
            R2 = np.zeros(len(T))
            for i, T_mid in enumerate(T):
                #Need at least 5 points to fit the polynomial
                if i > 5 and i < (len(T)-6):
                    #Separate the temperature and heat capacities into low and high range
                    (R2[i], a_low, a_high) = self._get_CpoR_R2(T, CpoR, i)
            max_R2 = max(R2)
            max_i = np.where(max_R2 == R2)[0][0]
            (max_R2, a_low_rev, a_high_rev) = self._get_CpoR_R2(T, CpoR, max_i)
            empty_arr = np.array([0.]*2)
            self.T_mid = T[max_i]
            self.a_low = np.concatenate((a_low_rev[::-1], empty_arr))
            self.a_high = np.concatenate((a_high_rev[::-1], empty_arr))

    def _get_CpoR_R2(self, T, CpoR, i_mid):
        """
        Calculates the R2 polynomial regression value.
        Parameters
        ----------
            T - (N,) ndarray
                Temperatures (K) to fit the polynomial
            CpoR - (N,) ndarray
                Dimensionless heat capacities that correspond to T array
            i_mid - int
                Index that splits T and CpoR arrays into a lower and higher range
        Returns
        -------
            R2 - float
                R2 value resulting from NASA polynomial fit to T and CpoR
            p_low - (5,) ndarray
                Polynomial corresponding to lower range of data
            p_high - (5,) ndarray
                Polynomial corresponding to high range of data
        """
        T_low = T[:i_mid]
        CpoR_low = CpoR[:i_mid]
        T_high = T[i_mid:]
        CpoR_high = CpoR[i_mid:]
        #Fit the polynomial
        p_low = np.polyfit(x = T_low, y = CpoR_low, deg = 4)
        p_high = np.polyfit(x = T_high, y = CpoR_high, deg = 4)

        #Find the R2
        CpoR_low_fit = np.polyval(p_low, T_low)
        CpoR_high_fit = np.polyval(p_high, T_high)
        CpoR_fit = np.concatenate((CpoR_low_fit, CpoR_high_fit))
        CpoR_mean = np.mean(CpoR)
        ss_reg = np.sum((CpoR_fit - CpoR_mean)**2)
        ss_tot = np.sum((CpoR - CpoR_mean)**2)
        R2 = ss_reg / ss_tot

        return (R2, p_low, p_high)

    def _fit_HoRT(self, HoRT0, T0 = c.T0('K')):
        """
        Calculates the a6 parameter for the NASA polynomial.
        Parameters
        ----------
            HoRT0 - float
                Dimensionless enthalpy at reference temperature
            T0 - float
                Reference temperature (K)
        """
        T_mid = self.T_mid
        a6_low = (HoRT0 - self._custom_HoRT(T0, self.a_low))*T0
        a6_high = (HoRT0 - self._custom_HoRT(T0, self.a_high))*T0

        #Correcting for offset
        H_low_last_T = self._custom_HoRT(T_mid, self.a_low) + a6_low/T_mid
        H_high_first_T = self._custom_HoRT(T_mid, self.a_high) + a6_high/T_mid
        H_offset = H_low_last_T - H_high_first_T

        self.a_low[5] = a6_low
        self.a_high[5] = T_mid * (a6_high/T_mid + H_offset)

    def _fit_SoR(self, SoR0, T0 = c.T0('K')):
        """
        Calculates the a7 parameter for the NASA polynomial.
        Parameters
        ----------
            SoR0 - float
                Dimensionless entropy at reference temperature
            T0 - float
                Reference temperature (K)
        """
        T_mid = self.T_mid
        a7_low = SoR0 - self._custom_SoR(T = T0, a = self.a_low)
        a7_high = SoR0 - self._custom_SoR(T = T0, a = self.a_high)
        #Correcting for offset
        S_low_last_T = self._custom_SoR(T_mid, self.a_low) + a7_low
        S_high_first_T = self._custom_SoR(T_mid, self.a_high) + a7_high
        S_offset = S_low_last_T - S_high_first_T

        self.a_low[6] = a7_low
        self.a_high[6] = a7_high + S_offset

    def _custom_HoRT(self, T, a):
        """
        Calculates the dimensionless enthalpy ignoring the contribution from the a6 parameter.
        Parameters
        ----------
            T - float
                Temperature
            a - (7,) ndarray
                NASA coefficients
        Returns
        -------
            HoRT - float
                Dimensionless enthalpy at T ignoring the a6 parameter
        """
        T_arr = np.array([1., T/2., T**2/3., T**3/4., T**4/5., 0., 0.])
        return np.dot(a, T_arr)

    def _custom_SoR(self, T, a):
        """
        Calculates the dimensionless entropy ignoring the contribution from the a7 parameter.
        Parameters
        ----------
            T - float
                Temperature
            a - (7,) ndarray
                NASA coefficients
        Returns
            SoR - float
                Dimensionless entropy at T ignoring the a7 parameter
        """
        T_arr = np.array([np.log(T), T, T**2/2., T**3/3., T**4/4., 0., 0.])
        return np.dot(a, T_arr)