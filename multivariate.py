import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import shapiro

plt.style.use('ggplot')

linkage_methods = ['ward', 'complete', 'average', 'single']
iterating_culture = ['corn1', 'corn2', 'soy']

nm_mun = []

for culture in self.inter_cul:
    plt.figure(figsize=(17, 8))
    plt.subplots_adjust(wspace=0.37, hspace=0.22)
    plt.suptitle('$1^{st}$ Harvest Maize', fontsize=15, y=0.95) if culture == 'corn1' else \
    plt.suptitle('$2^{nd}$ Harvest Maize', fontsize=15, y=0.95) if culture == 'corn2' else \
    plt.suptitle('$1^{st}$ Harvest Soybean', fontsize=15, y=0.95)
    
    gray_scale = ListedColormap(['#ffffff', '#D3D3D3', '#A9A9A9', '#696969', '#363636', '#000000'])
    classes = [ 2, 3, 4, 5, 6, 7]
    legend_elements = [Patch(facecolor=gray_scale(i), edgecolor='black', label=label) for i, label in enumerate(range(6))]
    vmin, vmax = 1, 7
    for i, year in enumerate(range(2010, 2020)):
        if culture == 'corn2':
            txt = pd.read_csv(os.path.join(scores_mpca, 'pca_scores_'+culture + '_' + str(year) + '.txt'), sep = ' ', index_col = [0],
                             encoding="latin-1")
            
            ano_safra = str(year)
            ano_prod  = str(year)[-2:]+'_'+str(year+1)[-2:]
            culture_index=culture + '_' + str(year) + '_' + str(year+1)
        
        else:
            txt = pd.read_csv(scores_mpca+'\\pca_scores_'+culture + '_' + str(year) + '_' + str(year+1) + '.txt', sep = ' ',
                              index_col = [0], encoding="latin-1")
            
            ano_safra = str(year) + '_' + str(year+1)
            ano_prod  = str(year)[-2:]+'_'+str(year+1)[-2:]
            culture_index=culture + '_' + str(year) + '_' + str(year+1)
        
        txt = txt.drop(columns= ['cd_mun', 'mesorregion', "x", "y"])
        cluster = txt.set_index('nm_mun')
        nm_mun.append([culture_index, len(cluster)])
        k_values = range(2, 10)
        cophenet_correlations = []
        silhouette_scores = []
        
        distances = pdist(cluster)

        for method in linkage_methods:
            Z = linkage(cluster, method=method)
            coph_corr, _ = cophenet(Z, distances)
            cophenet_correlations.append(coph_corr)

        most_method = [ method for method, score in zip(linkage_methods, cophenet_correlations) 
                       if score == max(cophenet_correlations) ]
        
        for k in k_values:
            model_ward = AgglomerativeClustering(n_clusters=k, affinity='euclidean', linkage=most_method[0])
            model_ward.fit(cluster)
            labels = model_ward.labels_
            score = silhouette_score(cluster, labels)
            silhouette_scores.append(score)

        valor_k = [ index + 2 for index, item in enumerate(silhouette_scores) if item == max(silhouette_scores) ] 
        c = AgglomerativeClustering(n_clusters=valor_k[0], affinity='euclidean', linkage=most_method[0])
        txt['class_agg'] = c.fit_predict(cluster) + 1
        all_data = pd.merge(parana, txt, how='left', left_on='nm_mun', right_on='nm_mun')
        all_data.loc[:, 'class_agg'] = all_data.loc[:, 'class_agg']#.fillna(-1)
        ax = plt.subplot(2, 5, i + 1)
                              
        ax.grid(False)
        ax.set_facecolor('white')
        plt.xlabel('lon', fontsize=12)
        plt.ylabel('lat', fontsize=12)
        if culture == 'corn2':
            ax.set_title('Harvest Year ' + ano_safra, fontsize=14)
        else:
            ax.set_title('Harvest Year ' + ano_prod.split('_')[0] + '/' + ano_prod.split('_')[1], fontsize=14)
        all_data.plot(column='class_agg', cmap=gray_scale, vmin=vmin, vmax=vmax, edgecolor='dimgray', linewidth = 0.3, ax=ax)
        ax.spines['top'].set_linewidth(0.9)
        ax.spines['right'].set_linewidth(0.9)
        ax.spines['bottom'].set_linewidth(0.9)
        ax.spines['left'].set_linewidth(0.9)
        ax.spines['top'].set_edgecolor('black')
        ax.spines['right'].set_edgecolor('black')
        ax.spines['bottom'].set_edgecolor('black')
        ax.spines['left'].set_edgecolor('black')
    ax.set_facecolor('white')
    frame = plt.legend(title = 'Clusters', handles=legend_elements, bbox_to_anchor = (-1.30, -0.22), ncol=9, fontsize=12).get_frame()
    
    ax.annotate('N', xy=(-1, -0.25), xytext=(-1, -0.42),
                arrowprops=dict(facecolor='black', width=5, headwidth=15),
                ha='center', va='center', fontsize=20,
                xycoords=ax.transAxes)
    frame.set_facecolor('white')
    plt.savefig(dir_cluster_hie + '\\english_' + culture + '_cluster_hierarquico_.png', bbox_inches='tight')
    plt.savefig(maps_save + '\\' + culture + '_cluster.jpeg', bbox_inches='tight', dpi=300)    