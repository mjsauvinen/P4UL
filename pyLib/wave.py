import numpy as np
#import gc  # garbage collector
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

class Wave:
    nWaves = int() # Counts the number of waves.
    writeHeader = True
    def __init__(self, time):
        Wave.nWaves += 1
        self.on = True
        self.time = time
        self.waveId = Wave.nWaves
        self.x = []; self.U = []
        self.Umean = None; self.Umax = None; self.Umin  = None
        self.center = None; self.width = None
    
    def __del__(self):
        # clearing the object.
        Wave.nWaves -= 1
        
        
    def newTime(self):
        Wave.writeHeader=True
    
    def add(self, Ui, xi):
        #print " adding id {}: U={} and x={} ".format(self.waveId,Ui,xi)
        self.U.append(Ui); self.x.append(xi) 
        
    def close(self):
        self.on = False    
        
    def analyze(self):
        if( len(self.U) > 0 ):
            Ua = np.array(self.U); xa = np.array(self.x)
            self.Umean = Ua.mean()
            self.Umax  = Ua.max(); self.Umin = Ua.min()
            self.center = xa.mean(); self.width = xa.max()-xa.min()
        else:
            self.Umean=self.Umax=self.Umin=0.
            self.center=self.width=0.
        
    def output(self, ofile=None):
        if( Wave.writeHeader ):
            header = "#{0:<8}\t{1:<8}\t{2:<8}\t{3:<8}\t{4:<8}\t{5:<4}\t{6:<5}"\
                .format( "Center","Umean","Umax","Umin","Width","Id","time")
            if( ofile ): 
                ofile.write(header+"\n")
            print header
            Wave.writeHeader = False
            
        outStr = " {0:<7f}\t{1:<7f}\t{2:<7f}\t{3:<7f}\t{4:<7f}\t{5:<4d}\t{6:<5g}"\
            .format(self.center,self.Umean,self.Umax, self.Umin,\
                self.width,self.waveId,self.time)
        if( ofile ): 
            ofile.write(outStr+"\n"); print outStr
        else:        
            print outStr



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

def waveInformation( time, Ui, xi, Ulimit, ofile=None ):
    waves = []; waves.append(Wave(time)); waves[0].newTime()
    j = 0

    for i in range(len(Ui)):
        #Check whether we're above the limit.
        if( Ui[i] > Ulimit ):
            # Then check whether the current wave is active.
            if( not waves[j].on ):
                waves.append(Wave(time))
                j += 1
            
            waves[j].add(Ui[i], xi[i])
        
        if( Ui[i] <= Ulimit ):
            if( waves[j].on and len(waves[j].U) > 0 ):
                waves[j].close() # Close the wave.

    for w in waves:
        w.analyze()
        w.output(ofile)
        del w

    del waves
    
    # Garbage collecting ... to avoid excessive memory consumption.
    #n = gc.collect(); print "{}".format(n)
    #print "gc.garbage = {}".format(gc.garbage)
    #gc.garbage[0].set_next(None); del gc.garbage[:]
    

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #


