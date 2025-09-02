from temperatureSolver.heat2d import Heat2d
import time


# define temperature solver
# height, width, time, nodes(mesh res), alpha

solver = Heat2d(time = 30, nodes=100)



solver.setBoundaryCondition('flux', 848000)

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