import numpy as np
from temperatureSolver.heat2d import Heat2d
import time

# Python wrapper for mui
# https://github.com/MxUI/MUI/tree/master/wrappers/Python
import mui4py


# define temperature solver
solver = Heat2d(1.0, 1.0, 1, 300, 1)

# Initial conditions
solver.initialiseTempField(0) # Set every cell to 20 degrees Celcius (293K)
solver.T[0, :] = 10 # Set left wall to 10 degrees

# MUI setup
appcomm = mui4py.mpi_split_by_app()
rank = appcomm.Get_rank()

dims = 2
config = mui4py.Config(dims, mui4py.FLOAT64)

domain = "tempSolver"
iface =  ["ifs1"]
unifaces = mui4py.create_unifaces(domain, iface, config)

uniface = unifaces["ifs1"]

uniface.set_data_types({"data": mui4py.FLOAT,
                       "temp": mui4py.FLOAT})

s_sampler = mui4py.SamplerExact()
t_sampler = mui4py.TemporalSamplerExact()


# Fetch initial conditions from ping
fetch_vals = uniface.fetch_many("data", solver.points, 0, s_sampler, t_sampler)
print(fetch_vals)

# Solver loop
i = 0

print("Starting timer, timestep : {:.2f}".format(solver.dt))
startTime = time.time()
while i < solver.time:

    solver.calculateHeatEquationWithNumpy(i)
    
    #uniface.push_many("temp", solver.points, solver.muiArr)
    i += solver.dt

runtime= time.time()-startTime
print("loop took {:.3f}".format(runtime))

solver.plotTemperature()