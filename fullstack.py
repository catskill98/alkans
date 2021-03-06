#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 21:38:43 2019

@author: antosha
"""

import openbabel as ob
#import warnings
#warnings.filterwarnings("ignore")
import pandas as pd  
import numpy as np  

import matplotlib.pyplot as plt  
import seaborn as seabornInstance 

from sklearn import metrics
from sklearn.model_selection import RepeatedKFold

from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import train_test_split 

from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Lasso
from sklearn.linear_model import LassoCV
from sklearn.linear_model import BayesianRidge
from sklearn.linear_model import Ridge
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV

from sklearn.decomposition import PCA

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

                i=0
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
###############################################################################
def fourregs(X, y, regressor, n):
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, 
                                                    test_size=0.2, 
                                                    random_state=0)
    regressor.fit(X_train, y_train)
#    y_pred = regressor.predict(X_test)
#    df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
    #df.plot(kind='bar',figsize=(10,8))
    #plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
    #plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
    #plt.show()
    print("R^2: ",regressor.score(X_test, y_test))
    
    
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, 
                                                        test_size=0.2, 
                                                        random_state=0)
    pca = PCA(n_components=n)
    pca.fit(X_train)
    X_train = pca.transform(X_train)
    X_test = pca.transform(X_test)
    regressor.fit(X_train, y_train)
#    y_pred = regressor.predict(X_test)
#    df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
    print ("Score with PCA: %.6f" % regressor.score(X_test, y_test))

def forregCV(X, y, regressor):
    kf = RepeatedKFold(n_splits=5, n_repeats=30)
   
    if regressor == ridge_regression:
        print("ridge")
        regressor = RidgeCV(alphas=[1e-3, 1e-2, 1e-1, 1, 10], cv=kf).fit(X, y)
        print('Ridge CV: ',regressor.score(X, y))
        return
    if regressor == lasso_regression:
        print("lasso")
        regressor = LassoCV(alphas=[1e-3, 1e-2, 1e-1, 1, 10], cv=kf).fit(X, y)
        print('Lasso CV: ',regressor.score(X, y))
        return
    #REPEATED KFOLD CROSS-VAL
#    midl = []
    scores = cross_val_score(regressor, X, y, cv=kf, n_jobs=-1)
    print("---> ", scores.mean())
#    for train_index, test_index in kf.split(X):
#    #    print("TRAIN:", train_index, "TEST:", test_index)
#        X_train, X_test = X[train_index], X[test_index]
#        y_train, y_test = y[train_index], y[test_index]
#        if regressor == ridge_regression:
#            regressor = RidgeCV(alphas=[1e-3, 1e-2, 1e-1, 1, 10]).fit(X, y)
#        elif regressor == lasso_regression:
#            regressor = LassoCV(alphas=[1e-3, 1e-2, 1e-1, 1, 10]).fit(X, y)
#        elif True:
#            regressor.fit(X_train, y_train)
##        y_pred = regressor.predict(X_test)
##        df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
#        midl.append(regressor.score(X_test, y_test))
#    #    print("linear splitted",regressor.score(X_test, y_test))
#    print('Average with repeatedKFOLD: ',  np.mean(midl))    
    
###############################################################################


print('Do you need BONDTYPE-marker?[0,1]')
marker = input()
print('input filename(only sdf format)')
na = input()
print('input parametr \n ALKAN: bp, \n GLASS: Tg (K)exp \n CANCERF4: A1')
desc = input()
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
a1 = []          #lists of fragments
a2 = []
a3 = []
etab = ob.OBElementTable();
i=0
n2=0
n3=0
while notatend:
    n2=(Alph2(obmol,f2,a2,marker,n2))
    n3=(Alph3(obmol,f3,a3,marker,n3))
    notatend = obconversion.Read(obmol)

notatend = obconversion.ReadFile(obmol,na)
i=0
n = 0
while notatend:
    Alph22(obmol,f2,a2,marker,n2)
    f2r.write(obmol.GetFormula()+"\n")
#    print(obmol.GetFormula()+"\n")
    notatend = obconversion.Read(obmol)
    n+=1
for element in a2:
    if len(etab.GetSymbol(element.name0)) == 1:
        element.name0 = etab.GetSymbol(element.name0)+'_'
    if len(etab.GetSymbol(element.name1)) == 1:
        element.name1 = etab.GetSymbol(element.name1)+'_'
    f2c.write(str(element.name0)+str(element.p0)+str(element.b0)+
             str(element.name1)+str(element.p1)+str(element.b1)+' ')
    
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
f2.close()
f3.close()
f2r.close()
f2c.close()
f3c.close()
f3r.close() 



abp = []
amp = []
i = 0


with open(na, 'r') as file:
    for line in file:
        if i == 1:
            abp.append(float(line.strip()))
            i = 0
            continue 
        if desc in line and i == 0:
            i = 1
#na = "pirimidines205_activity.txt"
#with open(na, 'r') as file:
#    for line in file:
#        abp.append(int(line.strip()))
           
       
ridge_regression = Ridge(alpha=0.1)  # alpha — величина регуляризации
regressor = LinearRegression()
lasso_regression = Lasso(alpha=0.1)
bayesian_regression = BayesianRidge()
n_pca = 5

midl = []
dataset = pd.read_csv('2matr.txt', sep = '\s+', header=None, engine = 'python')
print(dataset.shape)

from sklearn.metrics import silhouette_score



from sklearn.cluster import DBSCAN

dataset.describe()
X = dataset[:].values
y = pd.Series(abp)


X_pca = PCA(n_components=2).fit_transform(X)
plt.scatter(X_pca[:,0], X_pca[:,1])
dbscan = DBSCAN(eps=0.04, min_samples=2, metric="cosine")
res = dbscan.fit_predict(X_pca)
cluster_structure = res
plt.scatter(X_pca[:, 0], X_pca[:, 1], c=res, cmap=plt.cm.Paired)


print(silhouette_score(X,res))
#print(cluster_structure)
for i in np.unique(cluster_structure):
    if i == -1:
        continue
    X_cl = X[np.where(res==i)[0],:]
    y_cl = y[np.where(res==i)[0]]
    forregCV(X_cl, y_cl, regressor)
#    forregCV(X_cl, y_cl, ridge_regression)
#    forregCV(X_cl, y_cl, lasso_regression)
#    forregCV(X_cl, y_cl, bayesian_regression)
    
    
    
#    
#y_new = []
#ans_all = []
#iterat = 1
#for i in np.unique(cluster_structure):
#    if i == -1:
#        continue
#    X_in_clust = X[np.where(res==i)[0],:]
#    y_in_clust = y[np.where(res==i)[0],0]
#    ans = []
#    midl = []
#    for train_index, test_index in RepeatedKFold(min(10, X_in_clust.shape[0])).split(X_in_clust):
#        
#        
#        X_train = X_in_clust[train_index, :]
#        y_train = y_in_clust[train_index]
#
#        X_test = X_in_clust[test_index, :]
#        y_test = y_in_clust[test_index]
#        
#        
#        regressor.fit(X_train, y_train)
##        y_pred = regressor.predict(X_test)
##        df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
#        midl.append(regressor.score(X_test, y_test))
#        print("linear splitted",regressor.score(X_test, y_test))
#    

#from sklearn.neighbors import KNeighborsClassifier

#X_train, X_test, y_train, y_test = train_test_split(X, y, 
#                                                        test_size=0.2, 
#                                                        random_state=0)
##clustering = DBSCAN(eps=0.03, min_samples=2, n_jobs=-1)
#clustering = KNeighborsClassifier(n_neighbors=3)
##clustering.fit(X)
#clustering.fit(X_train, y_train)
#y_pred = clustering.predict(X_test)
#print("Accuracy:",metrics.accuracy_score(y_test, y_pred))
#print(res)
#print(clustering.score(X_train, y_train))
#print(res.labels_)
#print(silhouette_score(X_train,res))
#

#CLASSIC AND PCA 5
#for clf in [ridge_regression, regressor, lasso_regression, bayesian_regression]:
#    fourregs(X, y, clf, n_pca)
#    
##REPEATED KFOLD CROSS-VAL
#kf = RepeatedKFold(n_splits=5, n_repeats=5)
#for train_index, test_index in kf.split(X):
##    print("TRAIN:", train_index, "TEST:", test_index)
#    X_train, X_test = X[train_index], X[test_index]
#    y_train, y_test = y[train_index], y[test_index]
#    regressor.fit(X_train, y_train)
#    y_pred = regressor.predict(X_test)
#    df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
#    midl.append(regressor.score(X_test, y_test))
##    print("linear splitted",regressor.score(X_test, y_test))
#print('Average linear with KFOLD: ', np.mean(midl))    

#######################>>>trainCV regs<<<#######################
#for clf in [ridge_regression, regressor, bayesian_regression]:
#    forregCV(X, y, clf)


#######################>>>sh<<<#######################
#cv = ShuffleSplit(n_splits=5, test_size=0.2, random_state=0)
#print(cross_val_score(regressor, X, y, cv=5))
#cv = ShuffleSplit(n_splits=5, test_size=0.2, random_state=0)
#print(cross_val_score(ridge_regression, X, y, cv=5))
#cv = ShuffleSplit(n_splits=5, test_size=0.2, random_state=0)
#print(cross_val_score(lasso_regression, X, y, cv=5))
#cv = ShuffleSplit(n_splits=5, test_size=0.2, random_state=0)
#print(cross_val_score(bayesian_regression, X, y, cv=5))

#######################>>>PCA CV regs<<<#######################
#for clf in [regressor, ridge_regression,lasso_regression,bayesian_regression]:
#    for i in [1,2,3,4,5,6,7,8,9,10]:
#        
#        kf = RepeatedKFold(n_splits=5, n_repeats=5)
#        midl = []
#        for train_index, test_index in kf.split(X):
#            X_train, X_test = X[train_index], X[test_index]
#            y_train, y_test = y[train_index], y[test_index]
#            pca = PCA(n_components=i)
#            pca.fit(X_train)
#            X_train = pca.transform(X_train)
#            X_test = pca.transform(X_test)
#            clf.fit(X_train, y_train)
#            y_pred = clf.predict(X_test)
#            midl.append(clf.score(X_test, y_test))
#        print('Average with KFOLD i=',i,  np.mean(midl))
#    print("\n")
#    
       
#CROSS RIDGE AND LASSO
#clf = RidgeCV(alphas=[1e-3, 1e-2, 1e-1, 1, 10], cv=5).fit(X, y)
#print('Ridge CV: ',clf.score(X, y))
#clf = LassoCV(alphas=[1e-3, 1e-2, 1e-1, 1], cv=10).fit(X, y)
#print('Lasso CV: ',clf.score(X, y))


print('\n\n\n\n')
###############################################################################



dataset = pd.read_csv('3matr.txt', sep = '\s+', header=None, engine = 'python')
print(dataset.shape)
dataset.describe()
X = dataset[:].values
y = pd.Series(abp)



#X_pca = PCA(n_components=2).fit_transform(X)
#plt.scatter(X_pca[:,0], X_pca[:,1])
#clustering = DBSCAN(eps=0.3, min_samples=2, n_jobs=-1)
##clustering.fit(X)
#res = clustering.fit_predict(X)
#print(res)

#kf = RepeatedKFold(n_splits=5, n_repeats=5)
#midl = []
#for train_index, test_index in kf.split(X):
#    X_train, X_test = X[train_index], X[test_index]
#    y_train, y_test = y[train_index], y[test_index]
#    pca = PCA(n_components=2)
#    pca.fit(X_train)
#    X_train = pca.transform(X_train)
#    X_test = pca.transform(X_test)
##    clf.fit(X_train, y_train)
##    y_pred = clf.predict(X_test)
##    midl.append(clf.score(X_test, y_test))
#    dbscan = DBSCAN(eps=0.04, min_samples=2, metric="cosine")
#    res = dbscan.fit_predict(X_train)
#    cluster_structure = res
#    plt.scatter(X_train[:, 0], X_train[:, 1], c=res, cmap=plt.cm.Paired)
##print('Average with KFOLD i=',i,  np.mean(midl))
 
 
#print(res.labels_)
#print(silhouette_score(X,res))
 
 
 
#for clf in [regressor, bayesian_regression]:
#    for i in [1,2,3,4,5,6,7,8,9,10]:
#        
#        kf = RepeatedKFold(n_splits=5, n_repeats=20)
#        midl = []
#        for train_index, test_index in kf.split(X):
#            X_train, X_test = X[train_index], X[test_index]
#            y_train, y_test = y[train_index], y[test_index]
#            pca = PCA(n_components=i)
#            pca.fit(X_train)
#            X_train = pca.transform(X_train)
#            X_test = pca.transform(X_test)
#            clf.fit(X_train, y_train)
#            y_pred = clf.predict(X_test)
#            midl.append(clf.score(X_test, y_test))
#        print('Average with KFOLD i=',i,  np.mean(midl))
#    print("\n")

#for clf in [ridge_regression, regressor, lasso_regression, bayesian_regression]:
#    fourregs(X, y, clf, n_pca)
#
#
#kf = RepeatedKFold(n_splits=5, n_repeats=5)
#for train_index, test_index in kf.split(X):
##    print("TRAIN:", train_index, "TEST:", test_index)
#    X_train, X_test = X[train_index], X[test_index]
#    y_train, y_test = y[train_index], y[test_index]
#    regressor.fit(X_train, y_train)
#    y_pred = regressor.predict(X_test)
#    df = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
#    midl.append(regressor.score(X_test, y_test))
##    print("linear splitted",regressor.score(X_test, y_test))
#print('Average linear with KFOLD: ', np.mean(midl))  

#######################>>>trainCV regs<<<#######################
#for clf in [ridge_regression, regressor, bayesian_regression]:
#    forregCV(X, y, clf)

#cv = ShuffleSplit(n_splits=5, test_size=0.2, random_state=0)
#print(cross_val_score(regressor, X, y, cv=cv))
#cv = ShuffleSplit(n_splits=5, test_size=0.2, random_state=0)
#print(cross_val_score(ridge_regression, X, y, cv=cv))
#cv = ShuffleSplit(n_splits=5, test_size=0.2, random_state=0)
#print(cross_val_score(lasso_regression, X, y, cv=cv))
#cv = ShuffleSplit(n_splits=5, test_size=0.2, random_state=0)
#print(cross_val_score(bayesian_regression, X, y, cv=cv))
    
    
    
#for clf in [regressor, ridge_regression,lasso_regression,bayesian_regression]:
#    for i in [1,2,3,4,5,6,7,8,9,10]:
#        
#        kf = RepeatedKFold(n_splits=5, n_repeats=5)
#        midl = []
#        for train_index, test_index in kf.split(X):
#            X_train, X_test = X[train_index], X[test_index]
#            y_train, y_test = y[train_index], y[test_index]
#            pca = PCA(n_components=i)
#            pca.fit(X_train)
#            X_train = pca.transform(X_train)
#            X_test = pca.transform(X_test)
#            clf.fit(X_train, y_train)
#            y_pred = clf.predict(X_test)
#            midl.append(clf.score(X_test, y_test))
#        print('Average with KFOLD i=',i,  np.mean(midl))
#    print("\n")
#    
#    
#clf = RidgeCV(alphas=[1e-3, 1e-2, 1e-1, 1, 10], cv=5).fit(X, y)
#print('Ridge CV: ', clf.score(X, y))
#clf = LassoCV(alphas=[1e-3, 1e-2, 1e-1, 1], cv=10).fit(X, y)
#print('Lasso CV: ',clf.score(X, y))


del a2[:]
del a3[:]
