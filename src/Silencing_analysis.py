import math
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import norm, rankdata
from scipy.special import ndtr
from scipy.special import chdtrc as chi2_cdf
plt.rcParams.update({'font.size': 20})

color_list = ['green', 'orange', 'purple', 'gold', 'blue', 'red', 'brown']


#see https://en.wikipedia.org/wiki/Fisher%27s_method
#see https://en.wikipedia.org/wiki/Extensions_of_Fisher%27s_method

#implementation:
#https://academic.oup.com/bioinformatics/article/32/17/i430/2450768

#theory:
#https://doi.org/10.1016/S0167-7152(02)00310-3

#Brown's method for known covariance and Fishier's method

#for arbitrary dependancies, you can also check:
#https://www.nature.com/articles/s41598-022-07094-7


def main():
    analyse_inflammation()
    analyse_phenotypes()
    
    
def analyse_phenotypes():
    gene_knock = ['THBS1','FOSL1','CFH']
    phenotypes = ['Adhesion', 'Wound healing', 'Invasion']
    df = pd.read_excel('Data/Silencing_experiments/Silencing_data.xlsx', sheet_name='RAphenotypes')
    
    print(df)
    
    
    for gi,gene in enumerate(gene_knock):
           
        print()
        print()
        print(' Silencing Gene %s:' % gene)
        plt.figure(figsize=(6,2))
        plt.axhline(y=0, color='black', lw=1.5, ls='--')
        
        all_pval_fold = np.zeros((len(gene_knock),len(phenotypes)+1,3))
        for pi,pheno in enumerate(phenotypes):
            df_p = df[df['phenotype'] == pheno]
            
            fold = np.array(df_p['si%s' % gene] / df_p['siCTL'])
            fold = fold[~np.isnan(fold)]
            logfold = np.log2(fold)
            
            plt.scatter(pi*np.ones(len(logfold)), logfold, s=70, zorder=10, color=color_list[pi], edgecolors='black')
            plt.hlines(y=np.mean(logfold), xmin=pi-0.2, xmax=pi+0.2, linewidth=5, color=color_list[pi])
            
            all_pval_fold[gi,pi,:2] = stats.ttest_1samp(logfold, 0, alternative='less')
            all_pval_fold[gi,pi,2] = np.mean(logfold)

        all_pval_fold[gi,-1,0] = np.mean(all_pval_fold[gi,:-1,0])
        all_pval_fold[gi,-1,1] = combine_pvals_Fisher(all_pval_fold[gi,:-1,1])   
        all_pval_fold[gi,-1,2] = np.mean(all_pval_fold[gi,:-1,2])         
        plt.ylabel('log2fold')
        if gi==2:
            plt.xticks(np.arange(len(phenotypes)), phenotypes, rotation=90)
        else:
            plt.xticks([],[])
        plt.ylim(-3,3)
        #plt.xlabel('qPCR')
        plt.show()      
        
        print()        
        print('  folds: %s' % all_pval_fold[gi,:-1,2])
        print('  mean folds: %.2f'  % all_pval_fold[gi,-1,2])
        print()
        print('  tscores: %s' % all_pval_fold[gi,:-1,0])
        print('  mean tscore: %.2f'  % all_pval_fold[gi,-1,0])
        print()
        print('  p-vals: %s' % all_pval_fold[gi,:-1,1])
        print('  Combined (Fisher correction): %.6f' % all_pval_fold[gi,-1,1])
        
    


def analyse_inflammation():
    
    
    gene_to_study = ['IL6', 'CXCL10','MMP1','ACTA2','CDH11','C3']
    gene_knock = ['THBS1','FOSL1','CFH']
    covar_matrix = np.load('Data/Silencing_experiments/genes_covariance.npy')
    
    df = pd.read_excel('Data/Silencing_experiments/Silencing_data.xlsx', sheet_name='qPCR')
    df['fold'] = df['siRA'] / df['siCTL']
    df['log2fold'] = np.log2(df['fold'])
    
    all_pval_fold = np.zeros((len(gene_knock),len(gene_to_study)+1,3))
    for gi,gene in enumerate(gene_knock):
        
        print()
        print()
        print(' Silencing Gene %s:' % gene)
        plt.figure(figsize=(10,2))
        plt.axhline(y=0, color='black', lw=1.5, ls='--')
        
        df_g = df[df['gene_knocked'] == gene]
        for gj,gene_i in enumerate(gene_to_study):
            
            df_gi = df_g[df_g['gene_measured'] == gene_i]
    
            #gene_measured = df_g['gene_measured']
            #gene_knocked = df_g['gene_knocked']
            #cell_line = df_g['cell_line']
            siCTL = np.log2(df_gi['siCTL'].to_numpy())
            siRA = np.log2(df_gi['siRA'].to_numpy())
            logfold = df_gi['log2fold'].to_numpy()
            fold = df_gi['fold'].to_numpy()
            plt.scatter(gj*np.ones(len(logfold)), logfold, s=70, zorder=10, color=color_list[gj], edgecolors='black')
            plt.hlines(y=np.mean(logfold), xmin=gj-0.2, xmax=gj+0.2, linewidth=5, color=color_list[gj])
            
            #all_pval_fold[gi,gj,:2] = stats.wilcoxon(logfold, alternative='less')
            all_pval_fold[gi,gj,:2] = stats.ttest_rel(siRA, siCTL, alternative='less')
            #all_pval_fold[gi,gj,:2] = stats.ttest_1samp(logfold, 0, alternative='less')
            #all_pval_fold[gi,gj,:2] = stats.ttest_1samp(fold, 1, alternative='less')
            #This is a test for the null hypothesis that the expected value (mean) of a sample of independent observations is equal to the 0.
            all_pval_fold[gi,gj,2] = np.mean(logfold)
        
            
        all_pval_fold[gi,-1,0] = np.mean(all_pval_fold[gi,:-1,0])
        all_pval_fold[gi,-1,1] = combine_pvals_Brown(covar_matrix, all_pval_fold[gi,:-1,1])
        all_pval_fold[gi,-1,2] = np.mean(all_pval_fold[gi,:-1,2])         
        plt.ylabel('log2fold')
        if gi==2:
            plt.xticks(np.arange(len(gene_to_study)), gene_to_study, rotation=90)
        else:
            plt.xticks([],[])
        plt.ylim(-4,4)
        #plt.xlabel('qPCR')
        plt.show()
        print()        
        print('  folds: %s' % all_pval_fold[gi,:-1,2])
        print('  mean folds: %.2f'  % all_pval_fold[gi,-1,2])
        print()
        print('  tscores: %s' % all_pval_fold[gi,:-1,0])
        print('  mean tscore: %.2f'  % all_pval_fold[gi,-1,0])
        print()
        print('  p-vals: %s' % all_pval_fold[gi,:-1,1])
        print('  Combined (Brown correction): %.6f' % all_pval_fold[gi,-1,1])
        print('  Combined (independant): %.6f' % combine_pvals_Fisher(all_pval_fold[gi,:-1,1]))
        
        #print()
        #print('  p-vals: %s' % all_pval_fold[gi,:-1,1])
        #print('  Combined (Brown correction): %.4f' % all_pval_fold[gi,-1,1])
        #print('  Combined (independant): %.4f' % combine_pvals_Fisher(all_pval_fold[gi,:-1,1]))
        
    print()   
    print('Heatmaps')
    print(' t-scores')
    #plt.imshow(np.abs(all_pval_fold[:,:,0]), origin='lower', cmap='inferno', vmin=1)
    plt.imshow(-all_pval_fold[:,:,0], origin='lower', cmap='RdYlGn', vmin=0)
    plt.xticks(np.arange(len(gene_to_study)+1), gene_to_study + ['Comb.'], rotation=90)
    plt.yticks(np.arange(len(gene_knock)), gene_knock)
    plt.colorbar()
    plt.show()
    
    print()
    print(' p-vals')
    plt.imshow(all_pval_fold[:,:,1], origin='lower', cmap='inferno', vmin=0)
    #plt.imshow(all_pval_fold[:,:,1]<0.2, origin='lower', cmap='inferno')
    plt.xticks(np.arange(len(gene_to_study)+1), gene_to_study + ['Comb.'], rotation=90)
    plt.yticks(np.arange(len(gene_knock)), gene_knock)
    plt.colorbar()
    plt.show()
    
    
    sys.exit()
   
    
    p_values = [0.15, 0.15, 0.15, 0.15, 0.15, 0.15]
    covar_matrix = np.load('genes_covariance.npy')
    print(gene_to_study)
    print(p_values)
    #print(covar)
    print()
    print('Combined p-value Harmonic:', harmonic_mean(p_values)) #https://www.pnas.org/doi/10.1073/pnas.1814092116
    print("Combined p-value Fisher:", stats.combine_pvalues(p_values, method = 'fisher')[1])
    print("Combined p-value Stouffer:", stats.combine_pvalues(p_values, method = 'stouffer')[1])
    print()
    pval = combine_pvals_Brown(covar_matrix, p_values)
    print('Combined p-value Brown:', pval)
    pval_kov = combine_pvals_Brown(Kost_correlation(covar_matrix), p_values)
    print('Combined p-value Kost-Brown:', pval_kov)

    
    
def harmonic_mean(array, weights=None):
    
    if weights is None:
        weights = np.ones(len(array))/len(array)
        
    harm_mean = 1/(np.sum(weights/array))
    
    return harm_mean


def combine_pvals_Fisher(p_values):
    
    # the sum of k independent chi-squared values, each with two degrees of freedom,
    # follows a chi-squared distribution with 2k degrees of freedom.
    
    # Calculate the logarithm of each p-value
    log_p_values = np.log(p_values)
    
    # Sum the logarithms
    log_p_sum = np.sum(log_p_values)
    
    # Calculate the inverse chi-square distribution with degrees of freedom equal to twice the number of studies
    df = 2 * len(p_values)
    #chi2_inv = stats.chi2.ppf(q=0.95, df=df)
    
    # Multiply the result from step 2 by the result from step 3 to obtain the combined p-value
    p_combined = stats.chi2.sf(x=-2*log_p_sum, df=df)
    
    #print("P-values:", p_values)
    #print("Combined p-value:", p_combined)
    
    return p_combined
    
    
    
def combine_pvals_Brown(covar_matrix, p_values, extra_info = False):
    m = int(covar_matrix.shape[0])
    #print "m", m
    df_fisher = 2.0*m
    Expected = 2.0*m
    cov_sum = 0
    for i in range(m):
        for j in range(i+1, m):
            cov_sum += covar_matrix[i, j]
    
    #print "cov sum", cov_sum
    Var = 4.0*m+2*cov_sum
    c = Var/(2.0*Expected)
    df_brown = 2.0*Expected**2/Var
    if df_brown > df_fisher:
        df_brown = df_fisher
        c = 1.0

    x = 2.0*sum([-np.log(p) for p in p_values])
    #print "x", x
    p_brown = chi2_cdf(df_brown, 1.0*x/c)
    p_fisher = chi2_cdf(df_fisher, 1.0*x)
    
    if extra_info:
        return p_brown, p_fisher, c, df_brown
    else:
        return p_brown    
    
    
    
    
    
    
#Kost's approximation of the covariance between the -log cumulative distributions. This is calculated with a cubic polynomial fit.
def KostPolyFit(cor):
    a1, a2, a3 = 3.263, .710, .027 #Kost cubic coeficients
    return a1*cor+a2*cor**2+a3*cor**3

def Kost_correlation(cov_array):
    m = cov_array.shape[0]
    covar_matrix = np.zeros((m, m))
    for i in range(m):
        for j in range(i+1, m):
            cor = cov_array[i,j]
            covar = KostPolyFit(cor)
            covar_matrix[i, j] = covar
            covar_matrix[j, i] = covar
    return covar_matrix
    
    



if __name__ == '__main__':
    main()