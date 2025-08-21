import numpy as np
import matplotlib.pyplot as plt
import mui4py
from temperatureSolver.args import Args
from temperatureSolver.mui import MUI

class Heat2d:
    def __init__(self, time=1.0, dualSolver=False, maxTemp=100):
        # Parse solver arguments
        args = Args().args

        self.height = args.height # geometry
        self.width = args.length # geometry
        self.time = time # simulation time
        self.nodes = args.nodes # mesh res
        self.alpha = args.alpha # diffusivity
        self.solverNum = args.solverNum
        self.dualSolver = dualSolver
        print("Created heat solver with size {:f}m x {:f}m, {:d} nodes and diffusivity {:10f}".format(self.height, self.width, self.nodes**2, self.alpha))
        # derived properties
        self.dx = self.width/self.nodes
        self.dy = self.height/self.nodes    
         # if doing multiphase create a mui interface
        if self.dualSolver:
            self.numSolvers = 2
            self.mui = MUI(self.solverNum, self.nodes)
        else:
            self.numSolvers = 1
        # Define dt according to Courant-Friedrichs-Lewy condition in 2 dimensions
        # (https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition)
        self.dt = min(self.dx**2/(4*self.alpha), self.dy**2/(4*self.alpha))

        # if doing multisolver sync dt with the other interface
        if dualSolver:
            self.mui.uniface.push("dt", [self.solverNum, 0], self.dt)
            self.mui.uniface.push("alpha", [self.solverNum, 1], self.alpha)
            self.mui.commit( 0 )
            otherDt = self.mui.uniface.fetch("dt", [(self.solverNum+1)%2, 0], 0, mui4py.SamplerExact(), mui4py.TemporalSamplerExact())
            self.otherAlpha = self.mui.uniface.fetch("alpha", [(self.solverNum+1)%2, 1], 0, mui4py.SamplerExact(), mui4py.TemporalSamplerExact())
            self.alphaAvg = (self.otherAlpha + self.alpha)/(2*self.alpha)
            self.dt = min(self.dt, otherDt)
        else:
            self.otherAlpha = self.alpha
        print("dt set at: {:f}".format(self.dt))
        print("dx set at: {:f}".format(self.dx))
        print("dy set at: {:f}".format(self.dy))



        # Temperature array
        self.T = np.zeros((self.nodes, self.nodes))

        # Visualisation
        fig, self.axis = plt.subplots()
        self.pcm = self.axis.pcolormesh(self.T, cmap=plt.cm.jet, vmin=0, vmax = maxTemp)
        plt.colorbar(self.pcm, ax=self.axis)
        
        # setup for heat eq with numpy
        self.zerosX = np.zeros((1, self.nodes))
        self.zerosY = np.zeros((self.nodes, 1))
        self.zerosGhost = np.zeros((self.nodes, self.nodes-2))


    def initialiseTempField(self, val: float):
        self.T[:, :] = val

    def setTempIndex(self, x, y, val):
        if 0 <= x < self.nodes and 0 <= y < self.nodes:
            self.T[x][y] = val

    
    def calculateHeatEquation(self, time, animate = False):
        # local variable setup
        t = self.T.copy()
        xComponent = 0
        yComponent = 0
        fracX = self.dt/(self.dx**2)
        fracY = self.dt/(self.dy**2)

        # before doing anything fetch values from the left boundary (previous right boundary)
        # if self.solverNum == 0 there is no previous right boundary so ignore fetch and use initial condition
        if self.solverNum > 0:
            ## push T[0, :] (left boundary)
            self.mui.pushLeft(self.T[0, : ])
            self.mui.commit(time)
            #fetch right boundary of previous solver
            rightPrev = self.mui.fetchRightPrev(time)
            # set T[0, :] (left boundary ) according to forward difference approximation  
            # See README for more information on multilayer discretisation
            for y in range(0, self.nodes):
                xComponent = fracX * (self.alpha *(t[1][y]-t[0][y] + self.alphaAvg*(rightPrev[y]-t[0][y])))

                if  0 < y < self.nodes-1:
                    yComponent = fracY * (t[0][y+1] - 2*t[0][y] + t[0][y-1])
                # otherwise use a 'ghost flux' of 0
                elif y == 0:
                    yComponent = fracY * (t[0][y+1] - t[0][y])
                elif y == self.nodes-1:
                    yComponent = fracY * (t[0][y-1] - t[0][y])      
                 
                self.T[0][y] = t[0][y] + xComponent + self.alpha * yComponent
     
        
        for x in range(1, self.nodes):
            for y in range(0, self.nodes):
                # if we are not at the x boundaries use the regular forward differenece approximation for the x Component
                xComponent = None   
                if 0 < x < self.nodes-1:
                    xComponent = fracX * (t[x+1][y] - 2*t[x][y] + t[x-1][y])
                else:
                    #ghost flux
                    xComponent = fracX * (t[x-1][y] - t[x][y])
                
                # if we are not at the y boundaries use the regular forward differenece approximation for the y Component
                if  0 < y < self.nodes-1:
                    yComponent = fracY * (t[x][y+1] - 2*t[x][y] + t[x][y-1])
                # otherwise use a 'ghost flux' of 0
                elif y == 0:
                    yComponent = fracY * (t[x][y+1] - t[x][y])
                elif y == self.nodes-1:
                    yComponent = fracY * (t[x][y-1] - t[x][y])   
                        
                self.T[x][y] = t[x][y] + self.alpha * (xComponent + yComponent)

     
        #fetch left boundary of next solver if it exists
        leftNext = None
        if self.solverNum != self.numSolvers-1:
            leftNext = self.mui.fetchLeftNext(time)
        else:
            leftNext = np.zeros((self.nodes))
    
        # set T[-1, :] (left boundary ) according to forward difference approximation
        for y in range(0, self.nodes):
            xComponent = fracX * (self.alpha*(t[-2][y]-t[-1][y] + self.alphaAvg*(leftNext[y]-t[-1][y])))

            if  0 < y < self.nodes-1:
                yComponent = fracY * (t[-1][y+1] - 2*t[-1][y] + t[-1][y-1])
            # otherwise use a 'ghost flux' of 0
            elif y == 0:
                yComponent = fracY * (t[-1][y+1] - t[-1][y])
            elif y == self.nodes-1:
                yComponent = fracY * (t[-1][y-1] - t[-1][y])      
            
            self.T[-1][y] = t[-1][y] + xComponent +  self.alpha*yComponent
        
        # push the values at the rightBoundary
        if self.solverNum != self.numSolvers-1:
            self.mui.pushRight(self.T[-1, : ])
            self.mui.commit(time)
            
        if animate:
            self.pcm.set_array(self.T.T)
            self.axis.set_title("Temperature at time t: {:.3f}s".format(time))
            plt.pause(0.005)

    ## using numpy with shifts is roughly 150x faster for n=50
    ## use the above function to understand whats going on
    def calculateHeatEquationWithNumpy(self, time, animate = False):
        ## Idea: Do everything at once i.e xComponent =(dt/dx^2) Txplus1 + -2T + Txminus1
        rightPrev = None
        if self.solverNum > 0:
            ## push T[0, :] (left boundary)
            self.mui.pushLeft(self.T[0, : ])
            self.mui.commit(time)
            #fetch right boundary of previous solver
            rightPrev = [self.mui.fetchRightPrev(time)]
            rightPrev[0] *= self.alphaAvg

        leftNext = None
        if self.solverNum != self.numSolvers-1:
            leftNext = [self.mui.fetchLeftNext(time)]
            leftNext[0] *= self.alphaAvg

        Txminus1 = None
        Txplus1 = None
        Tx = self.T.copy()
        if self.solverNum == 0 and self.numSolvers == 1:
            Txplus1 = np.concatenate((self.T[:1], self.T[2:], self.T[-1:]))
            Txminus1 = np.concatenate((self.T[:1], self.T[:-2], self.T[-1:]))
        elif self.solverNum == 0:
            Txplus1 = np.concatenate((self.T[:1], self.T[2:], leftNext))
            Txminus1 = np.concatenate((self.T[:1], self.T[:-1]))
            Tx[-1, : ] *= self.alphaAvg
        elif self.solverNum == self.numSolvers - 1:
            Txplus1 = np.concatenate((self.T[1:], self.zerosX))
            Txminus1 = np.concatenate((rightPrev, self.T[:-1]))
            Tx[0, : ] *= self.alphaAvg
        else:
            #TODO: fix alpha stuff for more than 2 materials
            ## This is for if we have more than 2 solvers
            Txplus1 = np.concatenate((self.T[1:], leftNext))
            Txminus1 = np.concatenate((rightPrev, self.T[:-1]))
            Tx[-1, : ] *= self.alphaAvg
            Tx[0, : ] *= self.alphaAvg


        xComponent =  Txplus1 - Tx - self.T + Txminus1
        xComponent *= self.dt/(self.dx**2)
        # Y
        Typlus1 = np.concatenate((self.T[:, 1:], self.zerosY), axis=1)
        Tyminus1 = np.concatenate((self.zerosY, self.T[:, :-1]), axis=1)
        ghostFlux = np.concatenate((self.T[:, :1], self.zerosGhost, self.T[:, -1:]), axis=1)

        yComponent = Typlus1 -2*self.T + Tyminus1 + ghostFlux
        yComponent *= self.dt/(self.dx**2)
        self.T = self.T + self.alpha*(yComponent + xComponent)

        if self.solverNum != self.numSolvers-1:
            self.mui.pushRight(self.T[-1, : ])
            self.mui.commit(time)
        
        if animate:
            self.pcm.set_array(self.T.T)
            self.axis.set_title("Temperature at time t: {:.3f}s".format(time))
            plt.pause(0.005)

    def plotTemperature(self):

        self.pcm.set_array(self.T.T)
        self.axis.set_title("Temperature at time t: {:.3f}s".format(self.time))


        fig, ax = plt.subplots()

        xRange = np.arange(self.solverNum*self.width, (self.solverNum+1)*self.width, self.dx)
        tempData = np.average(self.T, axis=1)

        slope, intercept = np.polyfit(xRange, tempData, 1)

        ax.plot(xRange, tempData)
        ax.set_xlabel("x (m)")
        ax.set_ylabel("Temp (K)")
        ax.set_title("Temp from Solver {:d}. Equation of line: T = {:.3f}x + {:.2f}".format(self.solverNum, slope, intercept))
        
        if(self.solverNum == 0):
            print("T2 estimated as {:2f}".format(np.average(self.T[-1])))
        else:
            print("T2 estimated as {:2f}".format(np.average(self.T[0])))
            
        plt.show()

       