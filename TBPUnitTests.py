from TBPData import *
from TBPSimulate import *

def grabsolvers():
    solvers=[]
    for c in globals():
            c = globals()[c]
            if type(c)==type and len(c.mro())>2 and c.mro()[1]==Solver:
                solvers+=[c]
    return solvers

class Test:
    #Place stationary mass 1 particles at (-1,-1,-1) and (1,1,1) and a massless at (0,0,0)
    @staticmethod
    def acceltest1():
        ps = ParticleData()
        ps+=[1,[1,0,0],[0,0,0]]
        ps+=[5,[0,0,0],[0,0,0]]
        s = Solver(ps)
        acc = s.a(ps.state)
        expectacc = np.array([[-5*G,0,0],[G,0,0]])
        return np.isclose(acc,expectacc,atol=1e-15).all()

    @staticmethod
    def acceltest2():
        ps = ParticleData()
        ps+=[1,[-1,-1,-1],[0,0,0]]
        ps+=[1,[1,1,1],[0,0,0]]
        ps+=[0,[0,0,0],[0,0,0]]
        s = Solver(ps)
        acc = s.a(ps.state)
        expectacc = np.array([[G/12/np.sqrt(3)]*3,[-G/12/np.sqrt(3)]*3,[0]*3])
        return np.isclose(acc,expectacc,atol=1e-15).all()
    
    @staticmethod
    def zeroaccelevolve1():
        passed = True
        for dt in [1,0.1,1e-4,1e-8]:
            #Test all solvers to see if they handle it
            for c in grabsolvers():
                ps = ParticleData()
                ps+=[0,[0,0,0],[1,0,0]]
                ps+=[0,[4,5,6],[1,-5,20]]
                s = c(ps)
                tmp=s(dt)
                #Expected positions
                expect = np.array([[dt,0,0],[4+dt,5-5*dt,6+20*dt],[1,0,0],[1,-5,20]],dtype=np.float64)
                if not np.isclose(ps.state+tmp[1],expect,atol=1e-15).all():
                    passed=False
        return passed

    @staticmethod
    def testall():
        for f in Test.__dict__.keys():
            if type(Test.__dict__[f])==staticmethod and f!='testall':
                if not Test.__dict__[f]():
                    print(f"Test {f} failed")
                else:
                    print("Passed")

Test.testall()