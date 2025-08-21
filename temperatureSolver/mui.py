import mui4py
import numpy as np


class MUI:

    def __init__(self, solverNum, nodes):
        self.solverNum = solverNum
        self.nodes = nodes

        # Interface setup
        MUI_COMM_WORLD = mui4py.mpi_split_by_app()
        dims = 2
        config = mui4py.Config(dims, mui4py.FLOAT64)

        iface =  ["ifs"]
        domain = "Solver" + str(solverNum)
        unifaces = mui4py.create_unifaces(domain, iface, config)

        self.uniface = unifaces["ifs"]

        self.uniface.set_data_types({"temp": mui4py.FLOAT64,
                                     "dt": mui4py.FLOAT64,
                                     "alpha": mui4py.FLOAT64})

        #TODO: change spatial sampler so we can use different mesh resolution between solvers
        self.s_sampler = mui4py.SamplerExact()
        self.t_sampler = mui4py.TemporalSamplerExact()

        # rank = MUI_COMM_WORLD.Get_rank()
        # size = MUI_COMM_WORLD.Get_size()

        ## Point arrays for fetching using fetch_many
        self.fetchPointsLeft = np.zeros((self.nodes, 2))
        self.fetchPointsRight = np.zeros((self.nodes, 2))
        ## Point arrays for push using push_many
        self.pushPointsLeft = np.zeros((self.nodes, 2))
        self.pushPointsRight = np.zeros((self.nodes, 2))

        ## x-value within mui coord system of the boundary we want to fetch from (leftmost cell of solver to the right)
        fetchLeftBoundary = (self.nodes*(self.solverNum+1))
        #x-value within mui coord system of the boundary we want to push(leftmost cell of this solver)
        pushLeftBoundary = self.nodes*self.solverNum


        ## x-value within mui coord system of the boundary we want to fetch from (righmost most cell of solver to the left)
        fetchRightBoundary = (self.nodes*self.solverNum)-1
        #x-value within mui coord system of the boundary we want to push(rightmost cell of this solver)
        pushRightBoundary = self.nodes*(self.solverNum+1)-1

        for i in range(nodes):
            self.fetchPointsLeft[i] = [fetchLeftBoundary, i]
            self.fetchPointsRight[i] = [fetchRightBoundary, i]
            self.pushPointsLeft[i] = [pushLeftBoundary, i]
            self.pushPointsRight[i] = [pushRightBoundary, i]
    
    def pushRight(self, vals, data="temp"):
        self.uniface.push_many(data, self.pushPointsRight, vals)

    def pushLeft(self, vals, data="temp"):
        self.uniface.push_many(data, self.pushPointsLeft, vals)
    
    def commit(self, time):
        self.uniface.commit( time )

    def fetchRightPrev(self, time, data="temp"):
        vals = self.uniface.fetch_many(data, self.fetchPointsRight, time, 
                                       self.s_sampler, self.t_sampler)
        #print( vals )
        return vals
    
    def fetchLeftNext(self, time, data="temp"):
        vals = self.uniface.fetch_many(data, self.fetchPointsLeft, time, 
                                       self.s_sampler, self.t_sampler)
        #print( vals )
        return vals
