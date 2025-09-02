import numpy as np
import matplotlib.pyplot as plt
from temperatureSolver.args import Args
from temperatureSolver.mui import MUI

class Heat2d:
    def __init__(self, time=1.0, nodes=40):
        # Parse solver arguments
        args = Args().args

        self.height = args.height # geometry
        self.width = args.length # geometry
        self.time = time # simulation time
        self.nodes = nodes # mesh res
        self.alpha = args.alpha # diffusivity


        # derived properties
        self.dx = self.width/self.nodes
        self.dy = self.height/self.nodes 

        # Define dt according to Courant-Friedrichs-Lewy condition in 2 dimensions
        # (https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition)  
        self.dt = min(self.dx**2/(2*self.alpha), self.dy**2/(2*self.alpha))  

        # create MUI interface 
        self.mui = MUI( self.nodes, self.dt)

        # get solverNum and numSolvers from mui
        self.solverNum = self.mui.solverNum
        self.numSolvers = self.mui.numSolvers
        print("Solver: {:d}, Number of Solvers: {:d}".format(self.solverNum, self.numSolvers))

        print("Created heat solver {:d} with size {:f}m x {:f}m, {:d} nodes and diffusivity {:10f}".format(self.solverNum, self.height, self.width, self.nodes**2, self.alpha))

        
        # if doing multisolver create mui interface and sync dt and alpha with the other interfaces
        if self.numSolvers > 1:
            self.dt = self.mui.minDt
            
        prevAlpha, nextAlpha = self.mui.getAlphas(self.alpha)

        self.alphaAvgNext = (nextAlpha + self.alpha)/(2*self.alpha) 
        self.alphaAvgPrev = (prevAlpha + self.alpha)/(2*self.alpha)

        if self.solverNum == 0:
            print("dt set at: {:f}".format(self.dt))

        # Temperature array
        self.T = np.zeros((self.nodes, self.nodes))

        # Visualisation
        fig, self.axis = plt.subplots()
        self.pcm = self.axis.pcolormesh(self.T, cmap=plt.cm.jet, vmin=0, vmax = 100)
        plt.colorbar(self.pcm, ax=self.axis)

        # draw gridlines
        #self.axis.grid(which='major', axis='both', linestyle='-', color='0.8', linewidth=0.5)
        #self.axis.set_xticks(np.arange(0, self.nodes, 1))
        #self.axis.set_yticks(np.arange(0, self.nodes, 1))

        # setup for heat eq with numpy
        self.zerosX = np.zeros((1, self.nodes))
        self.zerosY = np.zeros((self.nodes, 1))
        self.zerosGhost = np.zeros((self.nodes, self.nodes-2))

        # set default BC
        if self.solverNum == 0:
            self.setBoundaryCondition('temp', 100)


    def initialiseTempField(self, val: float):
        self.T[:, :] = val

    def setBoundaryCondition(self, type: str, val: float, thermalConductivity=318):
        if type == 'flux':
            self.boundaryType = 'flux'
            self.boundaryValue = (np.full(shape=(1, self.nodes), fill_value=val, dtype=np.float64)*self.dx)/thermalConductivity
        elif type == 'temp':
            self.boundaryType = 'temp'
            self.boundaryValue = np.full(shape=(1, self.nodes), fill_value=val, dtype=np.float64)
            
        else:
            raise ValueError('Boundary condition can have type \'temp\' or \'flux\', not ' + type + ".")

    def setColorMapScale(self, min, max):
        self.pcm.set_clim(min, max)
       
        
    ## using numpy with shifts is roughly 150x faster for n=50

    def calculateHeatEquation(self, time, animate = False):
        ## Idea: Do everything at once i.e xComponent =(dt/dx^2) Txplus1 + -2T + Txminus1

        # push/fetch boundaries to/from adjacent solvers
        rightPrev = None
        if self.solverNum > 0:
            ## push T[0, :] (left boundary)
            self.mui.pushLeft(self.T[0, : ], time)
            #fetch right boundary of previous solver
            rightPrev = [self.mui.fetchRightPrev(time)]
            rightPrev[0] *= self.alphaAvgPrev
        else:
            if self.boundaryType == 'temp':
                rightPrev = self.boundaryValue
            else:
                rightPrev = self.boundaryValue + self.T[0, :]

        leftNext = None
        if self.solverNum != self.numSolvers-1:
            self.mui.pushRight(self.T[-1, : ], time)
            leftNext = [self.mui.fetchLeftNext(time)]
            leftNext[0] *= self.alphaAvgNext
        else:
            leftNext = self.zerosX

        
        Txminus1 = None
        Txplus1 = None
        Tx = self.T.copy()
        

        # set up shifted versions of T (xComponent)
        Txplus1 = np.concatenate((self.T[1:], leftNext))
        Txminus1 = np.concatenate((rightPrev, self.T[:-1]))

        # for edge solvers one of these is redundant, could be made more efficient
        if self.solverNum == 0:
            Tx[-1, : ] *= self.alphaAvgNext
        else:
            Tx[0, : ] *= self.alphaAvgPrev


    
        xComponent =  Txplus1 - Tx - self.T + Txminus1
        xComponent *= self.dt/(self.dx**2)

        # yComponent
        Typlus1 = np.concatenate((self.T[:, 1:], self.zerosY), axis=1)
        Tyminus1 = np.concatenate((self.zerosY, self.T[:, :-1]), axis=1)
        ghostFlux = np.concatenate((self.T[:, :1], self.zerosGhost, self.T[:, -1:]), axis=1)

        yComponent = Typlus1 -2*self.T + Tyminus1 + ghostFlux
        yComponent *= self.dt/(self.dy**2)
        self.T = self.T + self.alpha*(yComponent + xComponent)
        
        if animate:
            self.pcm.set_array(self.T.T)
            self.axis.set_title("Temperature at time t: {:.3f}s".format(time))
            plt.pause(0.005)

    def plotTemperature(self):

        self.pcm.set_array(self.T.T)
        self.axis.set_title("Temperature at time t: {:.3f}s".format(self.time))


        fig, ax = plt.subplots()

        xRange = np.linspace(self.solverNum*self.width, (self.solverNum+1)*self.width, self.nodes)
        tempData = np.average(self.T, axis=1)

        slope, intercept = np.polyfit(xRange, tempData, 1)

        ax.plot(xRange, tempData)
        ax.set_xlabel("x (m)")
        ax.set_ylabel("Temp (K)")
        ax.set_title("Temp from Solver {:d}. Equation of line: T = {:.3f}x + {:.2f}".format(self.solverNum, slope, intercept))
        if(self.solverNum != self.numSolvers - 1):
                print("T{:d} estimated as {:2f}".format(self.solverNum +2, np.average(self.T[-1])))
        if self.solverNum == 0 :
            if self.boundaryType == 'flux':
                print("T{:d} estimated as {:2f}".format(self.solverNum, np.average(self.T[0])))
            
        plt.show()

       