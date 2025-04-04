class StemGroup(object):
    """
       StemGroup Class
    """

    ########################################################################
    def __init__(self, extrj, maxDeviation=15, NStemlimit=25):

        # TODO
        self._PDBList = [];
        # TODO
        self._idStem = 0;
        # List of the number of stems per chain
        self._nStemsPerChainList = [];
        # List of the beta values after assigning values
        self._betaList = []
        # List of the occupation values after assigning values
        self._occList = []
        # TODO
        self._branchedStemsList = []
        # If   0 open and write headers for the very first time
        # else 1 open file to append
        self._writtenFLAG = 0
        # TODO
        self._maxDeviation = maxDeviation
        # TODO
        self._NStemlimit = NStemlimit
        # TODO
        self._nchains = 0
        # TODO
        self._exttrj = extrj
        # TODO
        self._numAtoms =  self._exttrj.num_atoms()
        # TODO
        self._grandlList = []
        # TODO
        self._grandnList = []
        # TODO
        self._grandsList = []
        # TODO
        self._frameList = []
        # TODO. These lists are only used in parallel processes
        self._ListStemAvg1 = []
        self._ListStemAvg2 = []
        self._ListStemAvg3 = []

    ########################################################################
    def look_for_stems(self, ini=0, end=-1, stride=1):

        """
        Look for stems in the frames from ini to end.

        Args:
            ini (int): First stem to look for. Default the fist one
            end (int): First stem to look for. Default the last one
            stride(int): Jump between frames.
                Example: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
                    ini=1, end=8, stride=2 => [1, 3, 5, 7]

        Returns:

        """

        # Steup the value of end by default
        if end == -1:
            end =  self._exttrj.get_numframes()

        for iframe in range(ini, end, stride):
            box = self._exttrj.universe.trajectory[iframe].dimensions[0:3]
            #print(self._exttrj.topology)
            #self._look_for_stem(iframe, box, isbranchArray, atbpArray)