import os
import sys
sys.path.insert(0,os.path.abspath('../../datatypes_and_database/'))
import admx_db_interface
from admx_db_interface import ADMXDB
from scipy.optimize import curve_fit
import numpy as np
import datetime
from datetime import datetime
from datetime import timedelta

def calc_sys_temp(scan_number, run="run1a"):
        """Function computes the system temperature. This function should be run on the fermilab computer, where the database is stored and where the ADMXDB class is defined."""

        #get digitizer scans, axion scans and NA logs
        if run=="run1a":
                port="5433"
                dbname="run1a"
        elif run=="run1b":
                port="5432"
                dbname="admx"
        print("Running on database {0} with port {1}".format(dbname, port))
        DB = ADMXDB(dbname=dbname, port=port)
        DB.hostname="admxdb01.fnal.gov"
        scan_number=str(scan_number)

        dig_table = "digitizer_log"
        what_to_select = "*"
        conditions=["digitizer_log_id=" + "'" + scan_number + "'"]
        dig_limit=10
        dig_query_string = DB.build_simple_select_query(from_table=dig_table, what_to_select=what_to_select, conditions=conditions, limit=dig_limit)
        dig_scan = DB.send_admxdb_query(query_string=dig_query_string)
        del dig_scan[0][7][0]

        from_table = "axion_scan_log"
        what_to_select="*"
        conditions=["digitizer_log_reference=" +"'"+ scan_number+"'"]
        limit=10
        query_string = DB.build_simple_select_query(from_table=from_table, what_to_select=what_to_select, conditions=conditions, limit=limit)
        axion_scan = DB.send_admxdb_query(query_string=query_string)

        naid = axion_scan[0][12]
        na_scan = DB.get_na_scan(naid = naid)

        #Get sensor values
        cavity_bottom_temp = axion_scan[0][8]
        cavity_top_temp = axion_scan[0][7]
        T_cavity = (cavity_bottom_temp + cavity_top_temp)/2
		
		condition1 = "sensor_name='msa_mk'"
		timestamp = dig_scan[0][0]
		timestamp_obj = datetime.strptime(timestamp, "%d-%m-%Y %I:%M:%S %p") 
		left_time_obj = timestamp_obj - timedelta(minutes=10)
		right_time_obj = timestamp_obj + timedelta(minutes=10)
		left_time = left_time_obj.strftime("%Y-%m-%d %I:%M:%S %p")
		right_time = right_time_obj.strftime("%Y-%m-%d %I:%M:%S %p")
		condition2 = "timestamp between "+"'"+left_time+"'"" and "+"'"+right_time+"'"
		
        query_string = DB.build_simple_select_query(from_table="sensors_log_double", what_to_select="*", conditions=[condition1, condition2], limit=10)
        temp_array = DB.send_admxdb_query(query_string=query_string) #THis is the temperature of A4 in the RF chain.
		time_array = [time.mktime(i[1].timetuple() for i in temp_array]
		value_array = [time.mktime(i[3].timetuple() for i in temp_array]
		squid_temperature = np.interp(time.mktime(timestamp_obj.timetuple()), time_array, value_array)
		
		fitted_func = power_fitter(dig_scan, axion_scan)
        R = np.max(fitted_func)/np.min(fitted_func)

        System_temp = (squid_temperature - T_cavity)/(R-1) #See supplementary material of the run1a results paper for equation explanation
        print("squid_temp"+str(squid_temperature), "T_cavity"+str(T_cavity))
        return {"System Temperature":System_temp, "Power Ratio":R}


def lorentzian(resonant_frequency, Q, size_of_lorentzian_in_bins, bin_size):
        length = size_of_lorentzian_in_bins
        bin = bin_size
        res_freq = resonant_frequency
        freqs = np.arange(start = res_freq - (length/2), stop = res_freq + (length/2), step = bin, dtype=np.dtype(float))

        lorentzian = (1/(2*np.pi))*(res_freq/Q)/((freqs-res_freq)**2 + (res_freq/(2*Q))**2)
        return lorentzian


def power_fitter(digitizer_scan, axion_scan):

        dig_scan = digitizer_scan
        power_spec = dig_scan[0][7] #power spectrum ch 1
        bin_size = dig_scan[0][5]*10**6   #frequency res ch1. in units of Hz
        Q = axion_scan[0][2] #Q ch1
        res_freq = axion_scan[0][3]*10**6     #mode freq ch1. in units of Hz
        length= len(power_spec)

        lorentzian = lambda f: (1/(2*np.pi))*(res_freq/Q)/((f-res_freq)**2 + (res_freq/(2*Q))**2)
        fitting_func = lambda f,a_0,a_1,a_2,a_3: (a_0 + a_1*lorentzian(f) + a_2*(f-res_freq)*lorentzian(f))*(1+a_3*(f-res_freq))

        xdata = np.arange(start = res_freq - (length/2)*bin_size, stop = res_freq + (length/2)*bin_size, step = bin_size, dtype=np.dtype(float))
        minimized_params, covariance = curve_fit(fitting_func, xdata, power_spec, p0=[np.max(power_spec),0,0,0], )
        return fitting_func(xdata, minimized_params[0], minimized_params[1], minimized_params[2], minimized_params[3])
