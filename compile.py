
import subprocess
from datetime import datetime
import os

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime(2009, 11, 11)
__modified__ = datetime(2012, 3, 14)
__version__  = "1.0"
__status__   = "Development, 11.11.2009, 14.3.2012"

def help():

    """
    @compile This is a simple script to call for automatic compiling of all
    fortran files necessary to run the soda2roms package. This is turned on in
    main.py with the compileAll=True

    Call this from command line using python compile.py
    """

def compileAll():
    logfile="compile.log"
    if os.path.exists(logfile): os.remove(logfile)
    log=open(logfile,'a')
    """Start the processes"""
    print "\n"

    print "Compiling IO.f90 to create ==> IOdata.so"
    proc = subprocess.Popen('f2py --verbose --fcompiler=intelem -c -m IOdata io.f90 --f90flags="-no-heap-arrays"',
                           shell=True, stdout=subprocess.PIPE,)
    stdout_value = proc.communicate()
    log.writelines(repr(stdout_value))

    print "Compiling behavior.f90 to create ==> behavior.so"
    proc = subprocess.Popen('f2py --verbose --fcompiler=intelem -c -m behavior behavior.f90 --f90flags="-no-heap-arrays"',
                           shell=True, stdout=subprocess.PIPE,)
    stdout_value = proc.communicate()[0]
    log.writelines(repr(stdout_value))

    print "Compiling calclight.f90 to create ==> calclight.so"
    proc = subprocess.Popen('f2py --verbose --fcompiler=intelem -c -m calclight calclight.f90 --f90flags="-no-heap-arrays"',
                           shell=True, stdout=subprocess.PIPE,)
    stdout_value = proc.communicate()[0]
    log.writelines(repr(stdout_value))

    print "Compiling perception.f90 to create ==> perception.so"
    proc = subprocess.Popen('f2py --verbose --fcompiler=intelem -c -m perception perception.f90 --f90flags="-no-heap-arrays"',
                           shell=True, stdout=subprocess.PIPE,)
    stdout_value = proc.communicate()[0]
    log.writelines(repr(stdout_value))

    print "Compiling predation.f90 to create ==> predation.so"
    proc = subprocess.Popen('f2py --verbose --fcompiler=intelem -c -m predF90 perception.f90 predation.f90 --f90flags="-no-heap-arrays"',
                           shell=True, stdout=subprocess.PIPE,)
    stdout_value = proc.communicate()[0]
    log.writelines(repr(stdout_value))

    print "Compiling zooplankton.f90 to create ==> zooplankton.so"
    proc = subprocess.Popen('f2py --verbose --fcompiler=intelem -c -m zooplankton zooplankton.f90 --f90flags="-no-heap-arrays"',
                           shell=True, stdout=subprocess.PIPE,)
    stdout_value = proc.communicate()[0]
    log.writelines(repr(stdout_value))

    print "Compiling bioenergetics.f90 to create ==> bioenergetics.so"
    proc = subprocess.Popen('f2py --verbose --fcompiler=intelem -c -m bioenergetics perception.f90 bioenergetics.f90 --f90flags="-no-heap-arrays"',
                           shell=True, stdout=subprocess.PIPE,)
    stdout_value = proc.communicate()[0]
    log.writelines(repr(stdout_value))

    log.close()

    print "Compilation finished and results written to file => %s"%(logfile)
    print "\n==================================================================="

compileAll()