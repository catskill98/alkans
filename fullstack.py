#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 21:38:43 2019

@author: antosha
"""
from sklearn.decomposition import PCA
import openbabel as ob
import pandas as pd  
import numpy as np  
import matplotlib.pyplot as plt  
import seaborn as seabornInstance 
from sklearn.model_selection import train_test_split 
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Lasso
from sklearn import metrics
from sklearn.linear_model import BayesianRidge
from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_score
'''
def printm(d,n,A,f):
    f.write("distance matr:\n")
    
    
    for i in range(n):
        #print(A[i].ind, A[i].name[:1], end = " ")  
        for j in range(n):
            f.write(str('{:>2}'.format(d[i][j]))+" ")
        f.write("\n") 
'''
class fragment(object):
    """class fragmentov"""
 
    def __init__(self, k, number, name0,p0,b0,name1,p1,b1,name2,p2,b2,
                 count
                 ):
        """Constructor"""
        self.k = k # 2 or 3 atoms 
        self.count = count # 
        self.num = number
        self.name0 = name0
        self.p0 = p0
        self.b0 = b0
        self.name1 = name1 
        self.p1 = p1
        self.b1 = b1
        self.name2 = name2 
        self.p2 = p2
        self.b2 = b2
###############################################################################
def getatom_b(obatom,marker):
    b=str(0)
    if (marker=="0"):
        return '*'
    if (obatom.HasSingleBond()):
        b='s'
    if (obatom.HasDoubleBond()):
        b='d'
    for obbond in ob.OBAtomBondIter(obatom):
        if (obbond.IsTriple()):
            b='t'
            break
    i=0
    for obbond in ob.OBAtomBondIter(obatom):
        if (obbond.IsDouble()):
            i+=1
    if (i>1):
        b='w'
    if (obatom.HasAromaticBond()):
        b='d'
        
    
    return b
###############################################################################
def Alph1(obmol,f,a1,marker,n):
    
    #print(n)
    #j=0
    
    for obatom in ob.OBMolAtomIter(obmol):
        i=1
        for element in a1:
            if (element.name0==obatom.GetAtomicNum() and
                element.p0==obatom.GetValence() and
                element.b0==getatom_b(obatom,marker)
                ):
                element.count+=1
                i=0
                break
        #print(i)       
        if (i==1):
            a1.append(fragment(1, n, obatom.GetAtomicNum(),obatom.GetValence(),
                               getatom_b(obatom,marker),0,0,0,  0,0,0,  1))
            #print(n)
            n+=1
            #print(n)
                    
    
    print("1-alph for mol ",obmol.GetFormula(), )
    for element in a1:
        print(element.num, "    ", element.name0, element.p0,element.b0,
              "   count =",element.count)
    return n
###############################################################################
def Alph11(obmol,f,a1,marker,n):
    n=0
    #print(n)
    #j=0
    
    for obatom in ob.OBMolAtomIter(obmol):
        i=1
        for element in a1:
            if (element.name0==obatom.GetAtomicNum() and
                element.p0==obatom.GetValence() and
                element.b0==getatom_b(obatom,marker)
                ):
                element.count+=1
                i=0
                break
        #print(i)       
        if (i==1):
            a1.append(fragment(1, n, obatom.GetAtomicNum(),obatom.GetValence(),
                               getatom_b(obatom,marker),0,0,0,  0,0,0,  1))
            #print(n)
            n+=1
            #print(n)
                    
    
    print("1-alph for mol ",obmol.GetFormula(), )
    for element in a1:
        print(element.num, "    ", element.name0, element.p0,element.b0,
              "   count =",element.count)
              
###############################################################################
def Alph2(obmol,f,a2,marker,n2):
    #n=0
    
    for obbond in ob.OBMolBondIter(obmol):
        i=1
        for element in a2:
            
            if ((element.name0==(obbond.GetBeginAtom()).GetAtomicNum()) and
            (element.p0==obbond.GetBeginAtom().GetValence()) and
            element.b0==getatom_b(obbond.GetBeginAtom(),marker)):
                
                
                if ((element.name1==(obbond.GetEndAtom()).GetAtomicNum()) and
                    (element.p1==obbond.GetEndAtom().GetValence()) and
                    element.b1==getatom_b(obbond.GetEndAtom(),marker)):
                    #element.count+=1
                    i=0
                    break
                
                
            if ((element.name0==(obbond.GetEndAtom()).GetAtomicNum()) and
            (element.p0==obbond.GetEndAtom().GetValence()) and
            element.b0==getatom_b(obbond.GetEndAtom(),marker)):
                
                
                if ((element.name1==(obbond.GetBeginAtom()).GetAtomicNum()) and
                    (element.p1==obbond.GetBeginAtom().GetValence()) and
                    element.b1==getatom_b(obbond.GetBeginAtom(),marker)):
                    #element.count+=1
                    i=0
                    break
            
        if (i==1):
            a2.append(fragment(2, n2, (obbond.GetBeginAtom()).GetAtomicNum(),
                               obbond.GetBeginAtom().GetValence(),
                               getatom_b(obbond.GetBeginAtom(),marker),
                               
                               (obbond.GetEndAtom()).GetAtomicNum(),
                               obbond.GetEndAtom().GetValence(),
                               getatom_b(obbond.GetEndAtom(),marker),
                               
                               0,0,0,1))
            n2+=1
    return n2
            
            
    
###############################################################################
def Alph22(obmol,f,a2,marker,n2):
    #n=0
    localmol=[]
    for i in range(n2):
        localmol.append(0)
    for obbond in ob.OBMolBondIter(obmol):
        i=1
        j=0
        for element in a2:
            
            if ((element.name0==(obbond.GetBeginAtom()).GetAtomicNum()) and
            (element.p0==obbond.GetBeginAtom().GetValence()) and
            element.b0==getatom_b(obbond.GetBeginAtom(),marker)):
                
                
                if ((element.name1==(obbond.GetEndAtom()).GetAtomicNum()) and
                    (element.p1==obbond.GetEndAtom().GetValence()) and
                    element.b1==getatom_b(obbond.GetEndAtom(),marker)):
                    #print('fo')
                    localmol[element.num]+=1
                    i=0
                    break
                
                
            if ((element.name0==(obbond.GetEndAtom()).GetAtomicNum()) and
            (element.p0==obbond.GetEndAtom().GetValence()) and
            element.b0==getatom_b(obbond.GetEndAtom(),marker)):
                
                
                if ((element.name1==(obbond.GetBeginAtom()).GetAtomicNum()) and
                    (element.p1==obbond.GetBeginAtom().GetValence()) and
                    element.b1==getatom_b(obbond.GetBeginAtom(),marker)):
                    #print('fo')
                    localmol[element.num]+=1
                    i=0
                    break
            j+=1
        if (i==1):
            pass
    #print(localmol)
    for i in range(n2):
        f.write(str(localmol[i])+" ")
    
    f.write("\n")
#        if (i==1):
#            a2.append(fragment(2, n, (obbond.GetBeginAtom()).GetType()[:1],
#                               obbond.GetBeginAtom().GetValence(),
#                               getatom_b(obbond.GetBeginAtom(),marker),
#                               
#                               (obbond.GetEndAtom()).GetType()[:1],
#                               obbond.GetEndAtom().GetValence(),
#                               getatom_b(obbond.GetEndAtom(),marker),
#                               
#                               0,0,0,1))
#            n+=1
            
            
            
#    print("2-alph for mol ",obmol.GetFormula())
#    for element in a2:
#        print(element.num, "  ", element.name0, element.p0,element.b0 ,
#              element.name1,element.p1,element.b1, " count =", element.count)
          
###############################################################################
def check1(a3,obmol,first,mutual,last,n3,localmol):
    i=1
    j=0
    for element in a3:
        if (element.name1==mutual.GetAtomicNum() and
            element.p1==mutual.GetValence() and
            element.b1==getatom_b(mutual,marker)):
            
            
            
            if (element.name0==first.GetAtomicNum() and
            element.p0==first.GetValence() and
            element.b0==getatom_b(first,marker) and
            
            
            element.name2==last.GetAtomicNum() and
            element.p2==last.GetValence() and
            element.b2==getatom_b(last,marker)):
           
           
                
               # element.name2==last.GetType()):
                
                
                i=0
                #localmol[j]+=1
                #element.count+=1
                break
                
            if (element.name2==first.GetAtomicNum() and
            element.p2==first.GetValence() and
            element.b2==getatom_b(first,marker) and
            
            
            element.name0==last.GetAtomicNum() and
            element.p0==last.GetValence() and
            element.b0==getatom_b(last,marker)):
           
                
                
                i=0
                #localmol[j]+=1
                #element.count+=1
                break
        j+=1
        
        
        
    if (i==1):
        a3.append(fragment(1, n3, first.GetAtomicNum(), first.GetValence(),
                               getatom_b(first,marker),
                               mutual.GetAtomicNum(), mutual.GetValence(),
                               getatom_b(mutual,marker),  last.GetAtomicNum(), 
                               last.GetValence(),
                               getatom_b(last,marker),  1))
        n3+=1
    return n3
        
        
###############################################################################
def check2(a3,obmol,first,mutual,last,n3,localmol):
    i=1
    j=0
    for element in a3:
        if (element.name1==mutual.GetAtomicNum() and
            element.p1==mutual.GetValence() and
            element.b1==getatom_b(mutual,marker)):
            
            
            
            if (element.name0==first.GetAtomicNum() and
            element.p0==first.GetValence() and
            element.b0==getatom_b(first,marker) and
            
            
            element.name2==last.GetAtomicNum() and
            element.p2==last.GetValence() and
            element.b2==getatom_b(last,marker)):
           
           
                
               # element.name2==last.GetType()):
                
                
                i=0
                localmol[element.num]+=1
                #element.count+=1
                
                
            if (element.name2==first.GetAtomicNum() and
            element.p2==first.GetValence() and
            element.b2==getatom_b(first,marker) and
            
            
            element.name0==last.GetAtomicNum() and
            element.p0==last.GetValence() and
            element.b0==getatom_b(last,marker)):
           
                
                
                i=0
                localmol[element.num]+=1
                #element.count+=1
                
        j+=1
        
        
        
    if (i==1):
        pass
    #print(localmol)
    return 
        
        
               
##ATOM --- GetId
##BOND --- GetIdx        
###############################################################################
def Alph3(obmol,f,a3,marker,n3):
    
    localmol=[]
    for i in range(n3):
        localmol.append(0)
   # n=0 
    calc=1
    for obbond in ob.OBMolBondIter(obmol):
        #for obbond in ob.OBMolBondIter(obmol):
        #print(obbond.GetIdx())
        for obbondin in ob.OBMolBondIter(obmol):
            if (obbondin.GetIdx()<calc):
                continue
            inend=(obbondin.GetEndAtom()).GetId()
            inbegin=(obbondin.GetBeginAtom()).GetId()
            end=(obbond.GetEndAtom()).GetId()
            begin=(obbond.GetBeginAtom()).GetId()
            #found joint atom 
            #1
            if (begin==inbegin):
                first=obbond.GetEndAtom()
                mutual=obbond.GetBeginAtom()
                last=obbondin.GetEndAtom()
                
                #   /list   /mol  /first -mutual- last
                n3=check1(a3,obmol,first,mutual,last,n3,localmol)
                continue
            #2
            if (begin==inend):
                first=obbond.GetEndAtom()
                mutual=obbond.GetBeginAtom()
                last=obbondin.GetBeginAtom()
                
                #   /list   /mol  /first -mutual- last
                n3=check1(a3,obmol,first,mutual,last,n3,localmol)
                continue
            #3
            if (end==inbegin):
                first=obbond.GetBeginAtom()
                mutual=obbond.GetEndAtom()
                last=obbondin.GetEndAtom()
                
                #   /list   /mol  /first -mutual- last
                n3=check1(a3,obmol,first,mutual,last,n3,localmol)
                continue
            #4
            if (end==inend):
                first=obbond.GetBeginAtom()
                mutual=obbond.GetEndAtom()
                last=obbondin.GetBeginAtom()
                
                #   /list   /mol  /first -mutual- last
                n3=check1(a3,obmol,first,mutual,last,n3,localmol)
                continue
            
        calc+=1
    return n3
###############################################################################
def Alph33(obmol,f,a3,marker,n3):
    
    localmol=[]
    for i in range(n3):
        localmol.append(0)
    calc=1
    for obbond in ob.OBMolBondIter(obmol):
        #for obbond in ob.OBMolBondIter(obmol):
        #print(obbond.GetIdx())
        for obbondin in ob.OBMolBondIter(obmol):
            if (obbondin.GetIdx()<calc):
                continue
            inend=(obbondin.GetEndAtom()).GetId()
            inbegin=(obbondin.GetBeginAtom()).GetId()
            end=(obbond.GetEndAtom()).GetId()
            begin=(obbond.GetBeginAtom()).GetId()
            #found joint atom 
            #1
            if (begin==inbegin):
                first=obbond.GetEndAtom()
                mutual=obbond.GetBeginAtom()
                last=obbondin.GetEndAtom()
                
                #   /list   /mol  /first -mutual- last
                check2(a3,obmol,first,mutual,last,n3,localmol)
                   
                continue
            #2
            if (begin==inend):
                first=obbond.GetEndAtom()
                mutual=obbond.GetBeginAtom()
                last=obbondin.GetBeginAtom()
                
                #   /list   /mol  /first -mutual- last
                check2(a3,obmol,first,mutual,last,n3,localmol)
                continue
            #3
            if (end==inbegin):
                first=obbond.GetBeginAtom()
                mutual=obbond.GetEndAtom()
                last=obbondin.GetEndAtom()
                
                #   /list   /mol  /first -mutual- last
                check2(a3,obmol,first,mutual,last,n3,localmol)
                continue
            #4
            if (end==inend):
                first=obbond.GetBeginAtom()
                mutual=obbond.GetEndAtom()
                last=obbondin.GetBeginAtom()
                
                #   /list   /mol  /first -mutual- last
                check2(a3,obmol,first,mutual,last,n3,localmol)
                continue
            
        calc+=1
    ##print(localmol)
    for i in range(n3):
        f.write(str(localmol[i])+" ")
    f.write("\n")
    #print("\n")
#    print("3-alph for mol ",obmol.GetFormula())
#    for element in a3:
#        print(element.num, "  ", 
#              element.name0,element.p0,element.b0,
#              element.name1,element.p1,element.b1,
#              element.name2,element.p2,element.b2," count =", element.count)
###############################################################################


print('Do you need BONDTYPE-marker?')
marker = input()
print('input filename(only sdf format')
na = input()
print('input parametr(like <bp>, <A1>, 3, 4')
desc = input()
#open file with mols
obconversion = ob.OBConversion()
obconversion.SetInFormat("sdf")
obmol = ob.OBMol()
f2=open('2matr.txt','w')
f3=open('3matr.txt','w')
f2r=open('2rows.txt','w')
f2c=open('2columns.txt','w')
f3r=open('3rows.txt','w')
f3c=open('3columns.txt','w')
notatend = obconversion.ReadFile(obmol,na)

#obconversion.SetInFormat("mol")
#notatend = obconversion.ReadFile(obmol,"2.mol")




a1 = []          #lists of fragments
a2 = []
a3 = []
etab = ob.OBElementTable();
#Alph1(obmol,f2,a1)
#mols in file
i=0
n2=0
n3=0
while notatend:
    #for obatom in ob.OBMolAtomIter(obmol):
    #    print(obatom.GetType()[:1], end=' ')
    #print("\n")
    #f.write('-------------MOL:   '+obmol.GetFormula()+'\n')
    
    
    #print(n)   CalcDesc(descnames=[])
    
    n2=(Alph2(obmol,f2,a2,marker,n2))
    n3=(Alph3(obmol,f3,a3,marker,n3))
    #print(obmol,n2, n3 )
    #Alph1(obmol,f2,a1,marker,0)
    #countfrag(obmol)
    #descriptmatrix
    notatend = obconversion.Read(obmol)
    #print("\n")
    #del a1[:]
    #del a2[:]
    #del a3[:]
    
    #a1.clear
    #a2.clear
    #a3.clear

#print(n2, n3, a2, a3)
#n2-=1
#n3+=1

notatend = obconversion.ReadFile(obmol,na)
i=0
n = 0
while notatend:
    #for obatom in ob.OBMolAtomIter(obmol):
    #    print(obatom.GetType()[:1], end=' ')
    #print("\n")
    #f.write('-------------MOL:   '+obmol.GetFormula()+'\n')
    
    #n=Alph1(obmol,f,a1,marker,n)-1
   
    
    Alph22(obmol,f2,a2,marker,n2)
    f2r.write(obmol.GetFormula()+"\n")
    #countfrag(obmol)
    #descriptmatrix
    notatend = obconversion.Read(obmol)
    #print("\n")
    #del a1[:]
    #del a2[:]
    #del a3[:]
    n+=1
    #print(n)
    #a1.clear
    #a2.clear
    #a3.clear
for element in a2:
    if len(etab.GetSymbol(element.name0)) == 1:
        element.name0 = etab.GetSymbol(element.name0)+'_'
    if len(etab.GetSymbol(element.name1)) == 1:
        element.name1 = etab.GetSymbol(element.name1)+'_'
    f2c.write(str(element.name0)+str(element.p0)+str(element.b0)+
             str(element.name1)+str(element.p1)+str(element.b1)+' ')

#    print(element.num, str(element.name0)+str(element.p0)+str(element.b0)+
#            str(element.name1)+str(element.p1)+str(element.b1)+
#            str(element.name2)+str(element.p2)+str(element.b2)+" ")

notatend = obconversion.ReadFile(obmol,na)
i=0
while notatend:
    f3r.write(obmol.GetFormula()+"\n")
    Alph33(obmol,f3,a3,marker,n3)
    notatend = obconversion.Read(obmol)
    
for element in a3:
    if len(etab.GetSymbol(element.name0)) == 1:
        element.name0 = etab.GetSymbol(element.name0)+'_'
    if len(etab.GetSymbol(element.name1)) == 1:
        element.name1 = etab.GetSymbol(element.name1)+'_'
    if len(etab.GetSymbol(element.name2)) == 1:
        element.name2 = etab.GetSymbol(element.name2)+'_'
        
    f3c.write(str(element.name0)+str(element.p0)+str(element.b0)+
            str(element.name1)+str(element.p1)+str(element.b1)+
            str(element.name2)+str(element.p2)+str(element.b2)+" ")
print("\n") 
#print("2-alph for mol ",obmol.GetFormula())

abp = []
amp = []
i = 0

#with open('ALKAN.SET','r') as file:
#    for line in file:
#        if 'A1' in line and i > 12:
#            abp.append(float((line.replace("  ", " ").split(' ')).pop(1)))
#        if 'A5' in line and i > 12:
#            amp.append(float((line.replace("  ", " ").split(' ')).pop(1)))
#           
#        i+=1
f2.close()
f3.close()
f2r.close()
f2c.close()
f3c.close()
f3r.close() 
with open(na, 'r') as file:
    for line in file:
        if i == 1:
            abp.append(float(line.strip()))
            i = 0
            continue 
        if desc in line and i == 0:
            i = 1
        
#if na == 'CANCERF4.SDF':
#    abp.pop[22]
#print(abp)
#print(len(abp), n, n2, n3)    




dataset = pd.read_csv('2matr.txt', sep = '\s+', header=None, engine = 'python')
print(dataset.shape)
dataset.describe()
X = dataset[:].values
#y = [[x] for x in abp]
y = pd.Series(abp)



X_train, X_test, y_train, y_test = train_test_split(X, y, 
                                                    test_size=0.2, 
                                                    random_state=0)
regressor = LinearRegression()  
regressor.fit(X_train, y_train)
y_pred = regressor.predict(X_test)
df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df.plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()
print("linear",regressor.score(X_test, y_test))



X_train, X_test, y_train, y_test = train_test_split(X, y, 
                                                    test_size=0.2, 
                                                    random_state=0)
pca = PCA(n_components=5)
pca.fit(X_train)
X_train = pca.transform(X_train)
X_test = pca.transform(X_test)
regressor.fit(X_train, y_train)
y_pred = regressor.predict(X_test)
df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df.plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()
print ("Score with PCA: %.6f" % regressor.score(X_test, y_test))


##############################

X_train, X_test, y_train, y_test = train_test_split(X, y, 
                                                    test_size=0.2, 
                                                    random_state=0)
ridge_regression = Ridge(alpha=0.1)  # alpha — величина регуляризации
ridge_regression.fit(X_train, y_train)
y_pred = ridge_regression.predict(X_test)
df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df.plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()
print('ridge  ', ridge_regression.score(X_test, y_test))

pca = PCA(n_components=5)
pca.fit(X_train)
X_train = pca.transform(X_train)
X_test = pca.transform(X_test)
ridge_regression.fit(X_train, y_train)
y_pred = ridge_regression.predict(X_test)
df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df.plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()
print ("Score with PCA: %.6f" % ridge_regression.score(X_test, y_test))

###################################
X_train, X_test, y_train, y_test = train_test_split(X, y, 
                                                    test_size=0.2, 
                                                    random_state=0)
lasso_regression = Lasso(alpha=0.1) # alpha — величина регуляризации
lasso_regression.fit(X_train, y_train)
y_pred = lasso_regression.predict(X_test)
df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df.plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()
print('lasso', lasso_regression.score(X_test, y_test))

pca = PCA(n_components=5)
pca.fit(X_train)
X_train = pca.transform(X_train)
X_test = pca.transform(X_test)
lasso_regression.fit(X_train, y_train)
y_pred = lasso_regression.predict(X_test)
df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df.plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()
print ("Score with PCA: %.6f" % lasso_regression.score(X_test, y_test))
###########################
X_train, X_test, y_train, y_test = train_test_split(X, y, 
                                                    test_size=0.2, 
                                                    random_state=0)
bayesian_regression = BayesianRidge()
bayesian_regression.fit(X_train, y_train)
y_pred = bayesian_regression.predict(X_test)
df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df.plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()
print('bayes', bayesian_regression.score(X_test, y_test))
pca = PCA(n_components=5)
pca.fit(X_train)
X_train = pca.transform(X_train)
X_test = pca.transform(X_test)
bayesian_regression.fit(X_train, y_train)
y_pred = bayesian_regression.predict(X_test)
df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df.plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()
print ("Score with PCA: %.6f" % bayesian_regression.score(X_test, y_test))
###########################
#plt.figure(figsize=(15,10))
#plt.tight_layout()
#seabornInstance.distplot(y[:])


#print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))  
#print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))  
#print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, y_pred)))


print('\n\n\n\n')
###############################################################################



dataset = pd.read_csv('3matr.txt', sep = '\s+', header=None, engine = 'python')
print(dataset.shape)
dataset.describe()
X = dataset[:].values
#y = [[x] for x in abp]
y = pd.Series(abp)
X_train, X_test, y_train, y_test = train_test_split(X, y, 
                                                    test_size=0.2, 
                                                    random_state=0)
regressor = LinearRegression()  
regressor.fit(X_train, y_train)
y_pred = regressor.predict(X_test)
df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df.plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()
print("linear",regressor.score(X_test, y_test))



X_train, X_test, y_train, y_test = train_test_split(X, y, 
                                                    test_size=0.2, 
                                                    random_state=0)
pca = PCA(n_components=5)
pca.fit(X_train)
X_train = pca.transform(X_train)
X_test = pca.transform(X_test)
regressor.fit(X_train, y_train)
y_pred = regressor.predict(X_test)
df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df.plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()
print ("Score with PCA: %.6f" % regressor.score(X_test, y_test))


##############################

X_train, X_test, y_train, y_test = train_test_split(X, y, 
                                                    test_size=0.2, 
                                                    random_state=0)
ridge_regression = Ridge(alpha=0.1)  # alpha — величина регуляризации
ridge_regression.fit(X_train, y_train)
y_pred = ridge_regression.predict(X_test)
df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df.plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()
print('ridge  ', ridge_regression.score(X_test, y_test))

pca = PCA(n_components=5)
pca.fit(X_train)
X_train = pca.transform(X_train)
X_test = pca.transform(X_test)
ridge_regression.fit(X_train, y_train)
y_pred = ridge_regression.predict(X_test)
df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df.plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()
print ("Score with PCA: %.6f" % ridge_regression.score(X_test, y_test))

###################################
X_train, X_test, y_train, y_test = train_test_split(X, y, 
                                                    test_size=0.2, 
                                                    random_state=0)
lasso_regression = Lasso(alpha=0.1) # alpha — величина регуляризации
lasso_regression.fit(X_train, y_train)
y_pred = lasso_regression.predict(X_test)
df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df.plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()
print('lasso', lasso_regression.score(X_test, y_test))

pca = PCA(n_components=5)
pca.fit(X_train)
X_train = pca.transform(X_train)
X_test = pca.transform(X_test)
lasso_regression.fit(X_train, y_train)
y_pred = lasso_regression.predict(X_test)
df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df.plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()
print ("Score with PCA: %.6f" % lasso_regression.score(X_test, y_test))
###########################
X_train, X_test, y_train, y_test = train_test_split(X, y, 
                                                    test_size=0.2, 
                                                    random_state=0)
bayesian_regression = BayesianRidge()
bayesian_regression.fit(X_train, y_train)
y_pred = bayesian_regression.predict(X_test)
df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df.plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()
print('bayes', bayesian_regression.score(X_test, y_test))
pca = PCA(n_components=5)
pca.fit(X_train)
X_train = pca.transform(X_train)
X_test = pca.transform(X_test)
bayesian_regression.fit(X_train, y_train)
y_pred = bayesian_regression.predict(X_test)
df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df.plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()
print ("Score with PCA: %.6f" % bayesian_regression.score(X_test, y_test))

#print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))  
#print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))  
#print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, y_pred)))








del a2[:]
del a3[:]
#fclose()
