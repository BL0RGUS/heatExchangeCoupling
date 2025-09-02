from temperatureSolver.heat2d import Heat2d
import time



# define temperature solver
# length, width, time, nodes(mesh res), alpha
solver = Heat2d(time=20.0, nodes=100)


# Solver loop
i = 0

print("Starting timer, timestep : {:.2f}".format(solver.dt))
startTime = time.time()
while i < solver.time:

    solver.calculateHeatEquationWithNumpy(i)
    i += solver.dt

runtime= time.time()-startTime
print("loop took {:.3f}".format(runtime))

solver.plotTemperature()