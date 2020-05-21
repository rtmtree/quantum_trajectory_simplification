import math
import random
import numpy as np
import time
import datetime
from projectq.ops import H, Z, X, Measure, All, Ph, Rz
from projectq.meta import Loop, Compute, Uncompute, Control


from projectq.backends import CircuitDrawer, CommandPrinter
from projectq.cengines import (MainEngine,
                               AutoReplacer,
                               LocalOptimizer,
                               TagRemover,
                               DecompositionRuleSet)
import projectq.setups.decompositions

# from hiq.projectq.cengines import GreedyScheduler, HiQMainEngine
# from hiq.projectq.backends import SimulatorMPI

from mpi4py import MPI
print('dd')
def createAllPossCaseNum(length,setBit):
    setPos=[]
    for i in range(0,pow(2,length)):
        if(count_set_bits(i)==setBit):
            setPos.append(Int2binString(i,length))
    return setPos
def createAllPossCaseSED(length):
    setPos=[]
    for i in range(0,pow(2,length)):
        setPos.append(Int2binString(i,length))
    return setPos
def cacl_SED(s,m,e):
    k = (m[2]-s[2])/(e[2]-s[2])
    px = ((1.0-k)*s[0]) + (k*e[0])
    py = ((1.0-k)*s[1]) + (k*e[1])
    sed = ((m[0]-px)*(m[0]-px))+((m[1]-py)*(m[1]-py))
    sqsed = np.sqrt(sed)
    return sqsed
def cacl_SUMSED(oriTrajectData,simTrajectData):
    sumSED = 0
    for i in range(1,len(simTrajectData)-1,1):
        if(simTrajectData[i]==True or simTrajectData[i]=="1"):
            # print("dont need to find SED")
            True
        else:
            sIndex =0
            eIndex =len(simTrajectData)-1
            for jj in range(i-1,0,-1):
                if(simTrajectData[jj]==True or simTrajectData[jj]=="1"):
                    sIndex=jj
                    break
            for j in range(i+1,len(simTrajectData),1):
                if(simTrajectData[j]==True or simTrajectData[j]=="1"):
                    eIndex=j
                    break
            s=oriTrajectData[sIndex]
            m=oriTrajectData[i]
            e=oriTrajectData[eIndex]
            sqsed = cacl_SED(s,m,e)
            # print(sIndex,i,eIndex,sqsed)
            sumSED+=sqsed
    return sumSED
def LSSEDSTATE(oriTrajectData,simTrajectData,memo=[]):
    newMemo=memo
    sumSED = 0
    if(len(oriTrajectData)!=len(simTrajectData)):
        print("oriTrajectData size dont match simTrajectData")
    for i in range(1,len(simTrajectData)-1,1):
        if(simTrajectData[i]==True or simTrajectData[i]=="1"):
            # print("dont need to find SED")
            True
        else:
            sIndex =i-1
            eIndex =i+1
            for jj in range(sIndex,-1,-1):
                if(simTrajectData[jj]==True or simTrajectData[jj]=="1"):
                    sIndex=jj
                    break
            for j in range(eIndex,len(simTrajectData),1):
                if(simTrajectData[j]==True or simTrajectData[j]=="1"):
                    eIndex=j
                    break
            # if
            # print(memo[sIndex])
            # if(memo!=[] and memo[sIndex][eIndex]!=False):
            #     print(memo[sIndex])
            #     print(sIndex,i,eIndex,'in memo')
            #     sqsed = memo[sIndex][eIndex]
            # else:
            s=oriTrajectData[sIndex]
            m=oriTrajectData[i]
            e=oriTrajectData[eIndex]
            sqsed = cacl_SED(s,m,e)
            # sqsed = cacl_SEDPAIR(oriTrajectData[sIndex:eIndex+1])
            if(memo!=[]):
                newMemo[sIndex][eIndex]=sqsed
            # print(sIndex,i,eIndex,sqsed)
            # print(sIndex,i,eIndex,sqsed)
            sumSED+=sqsed
    # print(newMemo)
    return [sumSED,newMemo]
def cacl_SEDPAIR(TrajectData):
    sumSED=0
    # print("1 to")
    # print(len(TrajectData)-2)
    for i in range(1,len(TrajectData)-1,1):
        # print("i")
        # print("come?")
        s=TrajectData[0]
        m=TrajectData[i]
        e=TrajectData[len(TrajectData)-1]
        sqsed = cacl_SED(s,m,e)
        sumSED+=sqsed
    return sumSED
def cacl_BESTNUM(oriTrajectData,setBit):
    bestPossCost = -1
    bestPossTra = [True] * len(oriTrajectData)
    for i in range(0,pow(2,len(oriTrajectData)-2),1):
        # print(count_set_bits(i))
        if(count_set_bits(i)==setBit-2):
            bools = Int2binString(i,len(oriTrajectData)-2)
            newBools = '1' + bools + '1'
            [curLSSED,memo]=LSSEDSTATE(oriTrajectData,newBools,[])
            # print(newBools,curLSSED)
            if(bestPossCost==-1 or curLSSED<bestPossCost):
                bestPossCost=curLSSED
                bestPossTra=newBools
    return [bestPossCost,bestPossTra]
def cacl_BESTSED(oriTrajectData,LSSEDthreshold):
    bestPossCost = 999999
    bestPossNum = len(oriTrajectData)
    bestPossTra = boolList2BinString([True] * len(oriTrajectData))
    for i in range(0,pow(2,len(oriTrajectData)-2),1):
        # print(count_set_bits(i))
        bools = Int2binString(i,len(oriTrajectData)-2)
        newBools = '1' + bools + '1'
        [curLSSED,memo]=LSSEDSTATE(oriTrajectData,newBools,[])
        # print(newBools , curLSSED)
        if(curLSSED<=LSSEDthreshold):
            # print(newBools , bestPossTra)
            # print(count_set_bitsString(newBools) , count_set_bitsString(bestPossTra))
            if(count_set_bitsString(bestPossTra)>count_set_bitsString(newBools) or (count_set_bitsString(bestPossTra)==count_set_bitsString(newBools) and bestPossCost>curLSSED)):
                bestPossCost=curLSSED
                bestPossTra=newBools
    print(bestPossCost, bestPossTra)
    return [bestPossCost,bestPossTra]
def intToBoolArray(integer,length):
    return [bool(integer & (1<<n)) for n in range(length)]

def count_set_bits(n):
    
    count = 0
    while n:
        n &= n - 1
        count += 1
    return count
def count_set_bitsString(n):
    count = 0
    if(n==0 or n==9999999):
        return 0
    for i in range(len(n)):
        if(n[i]=="1"):
            count+=1
    return count
def satisfyCheckFixN(binary,setbit):
    return count_set_bits(binString2Int(binary))==setbit
def satisfyCheckFixSED(binary,bestSED):
    return count_set_bitsString(binary)==count_set_bitsString(bestSED)
def toPercent(part,full):
    if(part.imag==0 and full.imag==0):
        # print("no imaginary part")
        return round((pow(part.real,2)/full.real) *100,10)
    return (pow(part,2)) *100
def boolList2BinString(lst):
    a = ''.join(['1' if x else '0' for x in lst])
    # return ''.join(reversed(a))
    return a
def binString2Int(lst):
    return int(lst, 2)
def Int2binString(lst,length):
    return boolList2BinString(intToBoolArray(lst,length))
def binString2BoolList(lst):
    return [True if n=="1" else False for n in lst]
def reverseBinary(bin):
    return ''.join(reversed(bin))
def grover_long(eng, Dataset, oracle, threshold,targetPoint=5,LSSEDthreshold=0,simType="#",prevSED=0,memo=[]):
    """
    Runs modified Grover's algorithm to find an element in a set larger than a 
    given value threshold, using the provided quantum oracle.

    Args:
        eng (MainEngine): Main compiler engine to run Grover on.
        Dataset(list): The set to search for an element.
        oracle (function): Function accepting the engine, an n-qubit register, 
        Dataset, the threshold, and an output qubit which is flipped by the
        oracle for the correct bit string.
        threshold: the threshold.

    Returns:
        solution (list<int>): the location of Solution.
    """
    N = pow(2,len(Dataset)-2)     
    n = int(math.ceil(math.log(N,2)))
    
    # number of iterations we have to run:
    num_it = int(math.sqrt(N)*9/4)
     
    m = 1
    landa = 6/5 
     
    # run num_it iterations
    i = 0
    while i<num_it:

        random.seed(i) 
        j = random.randint(0,int(m))
            
        x = eng.allocate_qureg(n)
    
        All(H) | x

        oracle_out = eng.allocate_qubit()
        X | oracle_out
        H | oracle_out
             
        with Loop(eng, j):

            oracle(eng, x, Dataset,threshold, oracle_out,targetPoint=targetPoint,simType=simType,LSSEDthreshold=LSSEDthreshold,prevSED=prevSED,memo=memo)
    
            with Compute(eng):
                All(H) | x
                All(X) | x
    
            with Control(eng, x[0:-1]):
                Z | x[-1]
    
            Uncompute(eng)
    
        All(Measure) | x
        Measure | oracle_out
         
        k=0
        xvalue=0
        z_i=""
        for i in range(len(x)):
            # print(int(x[k]))
            # xvalue=xvalue+int(x[k])*(2**k)
            if(int(x[i])==0):
                z_i+='0'
            else:
                z_i+='1'
            # k+=1
        # print(z_i)
        # print(int(x))
        # 
        #     print(int(x[i]))
        #     # if(int(x[i])==0):
        #     #     z_i+='0'
        #     # else:
        #     #     z_i+='1'
        [curLSSED,memo]=LSSEDSTATE(simTraj,'1'+z_i+'1',memo)
        if(simType=="#"):
            if(count_set_bitsString(z_i)==(targetPoint-2)):
                if( curLSSED<threshold ):
                    return [z_i,memo]
        else:
            if(((curLSSED<threshold and count_set_bitsString(z_i)==count_set_bitsString(prevSED)) or  count_set_bits(z_i)<count_set_bits(prevSED)) and curLSSED<LSSEDthreshold):
                return [z_i,memo]
            if(prevSED==0 and curLSSED<LSSEDthreshold):
                return [z_i,memo]
        m=m*landa
        i=i+1
        
    #fail to find a solution (with high probability the solution doesnot exist)
    return ["no answer",memo]
    
def oracleF( Dataset,d_0LSSED, simType="#",LSSEDthreshold=0,targetPoint=5,prevSED=0,memo=[]):
    # allPossibleCase=createAllPossCaseSED(len(Dataset)-2)
    # countMarked=0
    # fun = []
    # for elem in allPossibleCase:
    #     # print(elem)
    #     [curLSSED,memo] = LSSEDSTATE(simTraj,'1'+elem+'1',memo)
    #     if(simType=="#"):
    #         # print(elem)
    #         # print(count_set_bitsString(elem))
    #         if(count_set_bitsString(elem)==targetPoint and d_0LSSED>curLSSED):
    #             fun.append(1)
    #             countMarked+=1
    #         else:
    #             fun.append(0)
    #     else:
    #         if(count_set_bitsString(elem)<=count_set_bitsString(prevSED) and d_0LSSED>=curLSSED and LSSEDthreshold>=curLSSED):
    #             fun.append(1)
    #             countMarked+=1
    #         else:
    #             fun.append(0)
    if(simType=="#"):
        # print(elem)
        # print(count_set_bitsString(elem))
        def oracleMin(elem,d_0LSSED,memo):
            [curLSSED,memo] = LSSEDSTATE(simTraj,'1'+elem+'1',memo)
            
            if(count_set_bitsString(elem)==targetPoint and d_0LSSED>curLSSED):
                # print(elem,curLSSED,"mark")
                return 1
            else:
                # print(elem,curLSSED,"unmark")
                return 0
    else:
        # if(count_set_bitsString(elem)<=count_set_bitsString(prevSED) and d_0LSSED>=curLSSED and LSSEDthreshold>=curLSSED):
        #     fun.append(1)
        #     countMarked+=1
        # else:
        #     fun.append(0)
        def oracleMin(elem,d_0LSSED,prevSetbit,LSSEDthreshold,memo):
            [curLSSED,memo] = LSSEDSTATE(simTraj,'1'+elem+'1',memo)
            # print(count_set_bitsString(elem))
            # print(prevSetbit)
            # print(LSSEDthreshold)
            # print(curLSSED)
            #or (count_set_bitsString(elem)==prevSetbit and d_0LSSED>curLSSED)
            if((count_set_bitsString(elem)<prevSetbit or (count_set_bitsString(elem)==prevSetbit and d_0LSSED>curLSSED))  and LSSEDthreshold>curLSSED and (d_0LSSED==0 or d_0LSSED>curLSSED)):
                return 1
            else:
                return 0
    return oracleMin

def oracle(eng, qubits, Dataset, n0, output,simType="#",LSSEDthreshold=0,targetPoint=5,prevSED=0,memo=[]):
    """
    Marks the solutions by flipping the output qubit.

    Args:
        eng (MainEngine): Main compiler engine the algorithm is being run on.
        qubits (Qureg): n-qubit quantum register Grover search is run on.
        Dataset(list): The dataset.
        n0: The threshold.
        output (Qubit): Output qubit to flip in order to mark the solution..
    """
    # print(simType)
    
    allPossibleCase=createAllPossCaseSED(len(Dataset)-2)

    # print(n0)
    # print(allPossibleCase)
    # fun = [0]*len(allPossibleCase)
    fun = []
    # fun = [(LSSEDSTATE(simTraj,'1'+elem+'1',memo) - n0) for elem in allPossibleCase]
    for elem in allPossibleCase:
        [curLSSED,memo] = LSSEDSTATE(simTraj,'1'+elem+'1',memo)
        # print(curLSSED)
        if(simType=="#"):
            if(count_set_bits(binString2Int(elem))==targetPoint and n0>=curLSSED):
                fun.append(1)
            else:
                fun.append(0)
        else:
            if(count_set_bitsString(elem)>prevSED and n0>=curLSSED):
                fun.append(1)
            else:
                fun.append(0)
    # print(fun)
    # print(fun)
    a = sum(fun)
    n = int(math.ceil(math.log(pow(2,len(Dataset)-2),2)))
    # n = len(allPossibleCase)
    # print("gooo")
    # print(a)
    # print(fun)
    
    while a>0:
        num=[0]*n
        # print('come?')
        # print(fun.tolist())
        p = fun.index(1)
        # print(p)
        fun[p] = 0
        i=0
        a=a-1
        # print("P",p)
        while p/2 != 0:
            # print(i)
            # print(len(num))
            num[i] = p % 2
            p = p//2
            i = i+1       
        a1 = sum(num)  
        num1=num
        while a1>0:
            p = num1.index(1)
            a1 = a1-1
            num1[p]=0
            X | qubits[p]
             
        with Control(eng, qubits):
            X | output

        a1=sum(num)
        while a1>0:
            p = num.index(1)
            a1 = a1-1
            num[p]=0
            X | qubits[p]
        # print("sheet")
        # print(a)
        
        
def tracjectory_simplify(eng, Dataset,targetPoint=5,LSSEDthreshold=0,simType="#"):
    """
    The algorithm to find the maximum data in a dataset. The length of the input
    dataset should be 2^x, x is a positive integer.
    
    Args:
        eng (MainEngine): Main compiler engine the algorithm is being run on.
        Dataset(list): The dataset.
    """
    
    d_1=9999999
    
    N=len(Dataset)
    qubitLength=N-2
    #all possible solutions since we fixed 2 xs
    Nsolution=pow(2,N-2)

    memo= [[False for i in range(len(Dataset))] for j in range(len(Dataset))]
    
    if(simType=="#"):
        d_0=""
        for i in range(0,((len(Dataset)-2)-(targetPoint-2)) ):
            d_0+='0'
        for i in range(0,targetPoint-2):
            d_0+='1'
    else:     
        d_0=""
        for i in range(0,(len(Dataset)-2) ):
            d_0+='1'
    
    [d_0LSSED,memo] = LSSEDSTATE(simTraj,'1'+d_0+'1',memo)
    d_1LSSED=d_0LSSED
    # c=math.floor(1/2*pow(2,N-2))
    c=10
    i=0
    for i in range(0,c):
        #F = oracle
        if(i==0):
            print("initial guessing is ",'1'+d_0+'1')
            print("with SED",d_0LSSED)
        oracleFunction=oracleF(Dataset,d_0LSSED,simType=simType,LSSEDthreshold=LSSEDthreshold,targetPoint=targetPoint-2,memo=memo)
        # print(F)
        while (simType=="#" and ((d_1==9999999) or (d_1LSSED>d_0LSSED or count_set_bitsString(d_1)!=targetPoint-2))) or (simType=="SED" and ((d_1==9999999) or ((d_1LSSED>d_0LSSED and d_0LSSED!=0) or LSSEDthreshold<d_1LSSED or count_set_bitsString(d_1)> count_set_bitsString(d_0) ))):
            # print('d_1 > d_0')
            # print(simType)
            # print(d_0)
            # print(d_1)
            # print(d_1LSSED)
            # print(d_0LSSED)
            # print(count_set_bitsString(d_1))
            # print(count_set_bitsString(d_0))
            #map F to quantum state
            #grover long 
            # sineBeta = np.sqrt(Msolution/Nsolution)
            # beta = np.arcsin(sineBeta) 
            # # find J >= floor((np.pi-beta)/beta) + 1
            # J = math.floor((math.pi-beta)/beta) + 1
            # print("J",J)
            # angleRotation = 2* np.arcsin(  math.sin( np.pi/((4*J)+2) )  /  sineBeta  )
            # print("angleRotation",angleRotation)
            #phase rotation on qubits with J iteration
            if(simType=="#"):
                mapAllstate =[oracleFunction(Int2binString(elem,qubitLength),d_1LSSED,memo) for elem in range(pow(2,qubitLength))]
            else:
                mapAllstate =[oracleFunction(Int2binString(elem,qubitLength),d_1LSSED,count_set_bitsString(d_0),LSSEDthreshold,memo) for elem in range(pow(2,qubitLength))]
            # print(mapAllstate)
            if(sum(mapAllstate)==0):
                d_1=d_0
                break
            num=[0]*(qubitLength)
            p = mapAllstate.index(1)
            # print('p',p, LSSEDSTATE(Dataset,'1'+Int2binString(p,qubitLength)+'1') )
            mapAllstate[p] = 0
            mapAllstate[p] = 0
            binaryIndex=0
            while p/2 != 0:
                num[binaryIndex] = p % 2
                p = p//2
                binaryIndex+=1
            # print("num",num)
                # a1 = sum(elemnum)
                # num1=num.copy()
                # while a1>0:
                #     p = num1.index(1)
                #     a1 = a1-1
                #     num1[p]=0
                #     print("mark qubit at",p)
                #     X | z[p]
                # with Control(eng, z):
                #     X | oracle_out
                # print("num again",num)
                # a1=sum(num)
                # while a1>0:
                #     p = num.index(1)
                #     a1 = a1-1
                #     num[p]=0
                #     print("mark2 qubit at",p)
                #     X | z[p]
            z = eng.allocate_qureg(qubitLength)

            # start in uniform superposition
            All(H) | z

            # number of iterations we have to run:
            num_it = int(math.pi/4.*math.sqrt(1 << qubitLength))
            
            #phi is the parameter of modified oracle
            #varphi is the parameter of reflection across uniform superposition
            theta=math.asin(math.sqrt(1/(1 << qubitLength)))
            phi=math.acos(-math.cos(2*theta)/(math.sin(2*theta)*math.tan((2*num_it+1)*theta)))
            varphi=math.atan(1/(math.sin(2*theta)*math.sin(phi)*math.tan((2*num_it+1)*theta)))*2

            # prepare the oracle output qubit (the one that is flipped to indicate the
            # solution. start in state 1/sqrt(2) * (|0> - |1>) s.t. a bit-flip turns
            # into a (-1)-phase.
            oracle_out = eng.allocate_qubit()
            X | oracle_out
            H | oracle_out

            # run num_it iterations
            with Loop(eng, num_it):
                # oracle adds a (-1)-phase to the solution
                with Compute(eng):
                    #be 0
                    for k in range(qubitLength):
                        if(num[k]==0):
                            X | z[k]
                    # X | z[1]
                    # X | z[4]
                    # X | z[5]
                with Control(eng, z):
                    X | oracle_out
                Uncompute(eng)

                # reflection across uniform superposition
                with Compute(eng):
                    All(H) | z
                    All(X) | z

                with Control(eng, z[0:-1]):
                    Z | z[-1]

                Uncompute(eng)
        
            # prepare the oracle output qubit (the one that is flipped to indicate the
            # solution. start in state |1> s.t. a bit-flip turns into a e^(i*phi)-phase. 
            H | oracle_out
            with Compute(eng):
                #be 1
                for k in range(qubitLength):
                    if(num[k]==1):
                        X | z[k]
                # X| z[0]
                # X| z[2]
                # X| z[3]
            with Control(eng, z):
                Rz(phi) | oracle_out
                Ph(phi/2) | oracle_out

            Uncompute(eng)

            with Compute(eng):
                All(H) | z
                All(X) | z

            with Control(eng, z[0:-1]):
                Rz(varphi) | z[-1]
                Ph(varphi/2) | z[-1]

            Uncompute(eng)
                        

            All(Measure) | z
            Measure | oracle_out
            eng.flush()
            xprime_i=""
            for z_i in z:
                if(int(z_i)==0):
                    xprime_i=xprime_i+'0'
                else:
                    xprime_i=xprime_i+'1'
            xprime=xprime_i
            d_1= xprime
            print("d_1",d_1)
            [d_1LSSED,memo] = LSSEDSTATE(simTraj,'1'+d_1+'1',memo)
        if (simType=="#" and (d_1LSSED<d_0LSSED and count_set_bitsString(d_1)==targetPoint-2)) or (simType=="SED" and (d_1LSSED<d_0LSSED and count_set_bitsString(d_1)<count_set_bitsString(d_0))):
            print('d_1 < d_0')
            i=0
        elif (simType=="#" and (d_1LSSED==d_0LSSED and count_set_bitsString(d_1)==targetPoint-2)) or (simType=="SED" and (d_1LSSED<d_0LSSED and count_set_bitsString(d_1)==count_set_bitsString(d_0))):
            print('d_1 == d_0')
            print("d_1",d_1)
            # return '1'+d_0+'1'
            shit=1
            return '1'+d_1+'1'
        print('set d_0=d_1')
        d_0=d_1
        d_0LSSED=d_1LSSED
        if(simType=="#"):
            d_1=9999999
        else:
            d_1=9999999
    print("break with c",c)
    return '1'+d_0+'1'
    # [a,memo]=grover_long(eng, Dataset, oracle, GuesspathThreshold,targetPoint=targetPoint,LSSEDthreshold=LSSEDthreshold,simType=simType,prevSED=a1,memo=memo)
    # i=0
    # while i<int(math.log(pow(2,len(Dataset)-2), 2)):
    #     [GuesspathThreshold,memo] = LSSEDSTATE(simTraj,'1'+a1+'1',memo) 
    #     if a!="no answer":
    #         a1=a
    #         [a,memo]=grover_long(eng, Dataset, oracle, GuesspathThreshold,targetPoint=targetPoint,LSSEDthreshold=LSSEDthreshold,simType=simType,prevSED=a1,memo=memo)
    #         i=i+1
    #     else:
    #         [a,memo]=grover_long(eng, Dataset, oracle, GuesspathThreshold,targetPoint=targetPoint,LSSEDthreshold=LSSEDthreshold,simType=simType,prevSED=a1,memo=memo)
    #         i=i+0.5
    # return('1'+Int2binString(a1,len(Dataset)-2)+'1')

def _is_power2(num):
    """
    Checks whether a given number is a power of two
    """
    return num != 0 and ((num & (num - 1)) == 0)

if __name__ == "__main__":

    trajSize = 10
    targetPoint = 5 # add 2 for head and tail
    # targetPoint = math.floor(trajSize/10) # add 2 for head and tail
    LSSEDthreshold = 0.000125
    startFrom=10
    simType="#"
    fileName = "Geolife000.plt"

    f=open(fileName, "r").read()
    allTraj=f.split('\n')
    simTrajRaw =allTraj[startFrom:startFrom+trajSize]
    simTraj = []
    for i in range(0,len(simTrajRaw),1):
        eachRow = simTrajRaw[i].split(",")
        dateTime = eachRow[5]+" "+eachRow[6]
        # print(dateTime)
        timeStamp = time.mktime(datetime.datetime.strptime(dateTime, "%Y-%m-%d %H:%M:%S").timetuple())
        x = np.array([eachRow[0],eachRow[1],timeStamp])
        # print(x)
        simTraj.append(x.astype(np.float))
    countRound=0
    countSatisfyfixN=0
    countBestfixN=0
    countSatisfyfixSED=0
    countBestfixSED=0
    roundRepeat=100
    backend = SimulatorMPI(gate_fusion=True, num_local_qubits=20)
    cache_depth = 10
    rule_set = DecompositionRuleSet(modules=[projectq.setups.decompositions])
    engines = [TagRemover(),
    		   LocalOptimizer(cache_depth),
    		   AutoReplacer(rule_set),
    		   TagRemover(),
    		   LocalOptimizer(cache_depth),
    		   #,CommandPrinter(),
    		   GreedyScheduler()
    		   ]
    eng = HiQMainEngine(backend, engines) 

    Dataset = simTraj
    
    

    # if not _is_power2(len(Dataset)):
    #     raise ValueError("The size of the dataset must be a power of 2!")

    while(True):
        soluHighestProb = tracjectory_simplify(eng, Dataset,targetPoint,LSSEDthreshold=LSSEDthreshold,simType=simType)
        countRound=countRound+1
        print("=======================================================================")
        print("---------------------------------------RESULT---------------------------------------------")
        print("solution with highest prob is ", ((soluHighestProb)))
        [soluSED,memoBuff]=LSSEDSTATE(Dataset,soluHighestProb,[]) 
        print("SED of highest prob is ",soluSED)

        if(simType=="#"):
            bestSED,bestSol = cacl_BESTNUM(simTraj,targetPoint)
            print("Fixed # ---------------",targetPoint,"--------------")
            print("solution of best possible is ",(bestSol))
            print("SED of best possible is ", bestSED)
            if(satisfyCheckFixN((soluHighestProb),targetPoint)):
                countSatisfyfixN=countSatisfyfixN+1
                print("!!!!Congratulation!!!  the highest prob satisfy the requirement!!!!!!")
                if(bestSED==soluSED):
                    countBestfixN=countBestfixN+1
                    print("!!!!Congratulation!!!  you got the best one!!!!!!!!!")
            print ("Satisfactorty Rate fix# = ",countBestfixN,":",countSatisfyfixN,"/",countRound)
        else:
            bestSED,bestSol = cacl_BESTSED(simTraj,LSSEDthreshold)
            print("Fixed SED ---------------",LSSEDthreshold,"--------------")
            print("solution of best possible is ",bestSol)
            print("SED of best possible is ", bestSED)
            if(soluSED<LSSEDthreshold and satisfyCheckFixSED((soluHighestProb),(bestSol))):
                countSatisfyfixSED=countSatisfyfixSED+1
                print("!!!!Congratulation!!!  the highest prob satisfy the requirement!!!!!!")
                if(bestSED==soluSED):
                    countBestfixSED=countBestfixSED+1
                    print("!!!!Congratulation!!!  you got the best one!!!!!!!!!")
            print ("Satisfactorty Rate fixSED = ",countBestfixSED,":",countSatisfyfixSED,"/",countRound)


        if(countRound==roundRepeat and simType=="#"):
            simType="SED"
            countRound=0
            # break
        if(countRound==roundRepeat and simType=="SED"):
            print ("Satisfactorty Rate fix# = ",countBestfixN,":",countSatisfyfixN,"/",countRound)
            print ("Satisfactorty Rate fixSED = ",countBestfixSED,":",countSatisfyfixSED,"/",countRound)
            break