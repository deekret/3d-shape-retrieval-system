
import matplotlib.patches as mpatches
import matplotlib
import numpy as np
from sklearn import manifold
import seaborn as sns
import matplotlib.colors as mpl_colors
import matplotlib.pyplot as plt
import pandas as pd
import os

#Preprocessing data

features = pd.read_csv("outputStand.csv", skiprows=3,sep=";",index_col=False)
features[['name','mesh_id']] = features['name'].astype(str).str.split("/",expand=True)
features['mesh_id'] = features['mesh_id'].str.replace('.off', '')

# get_labels Returns a dictionary with class as key and mesh as value

def get_labels():
    """
    Generates dictionary with class label for each mesh,
    and returns it.
    """
    file_path = 'DB' #insert here path to your database
    classes = {}
    values = []
    for root, dirs, files in os.walk(file_path, topdown=True):
        for file in files:
            if file.endswith('.off'):
                index = file.split('.', 1)[0].replace('m', '')
                classes[index] = os.path.basename(root)
                
    return classes


def update_annot(scatter, annot, ind, norm, cmap, labels, mesh_ids, classes):
    """
    Update the tooltip with the correct mesh_id and class name.
    """
    idx = ind["ind"][0]
    mesh_id = mesh_ids[idx]

    pos = scatter.get_offsets()[idx]
    annot.xy = pos
    text = f"mesh: {mesh_id}, class: {classes[str(mesh_id)]}"
    annot.set_text(text)
    #annot.get_bbox_patch().set_facecolor(cmap(norm(labels[(row.iloc[0])])))
    annot.get_bbox_patch().set_alpha(0.6)


def hover(fig, scatter, annot, ax, event, norm, cmap,
          labels, mesh_ids, classes):
    """
    Implements a hover event that call the update annot function.
    Continuously checks whether the mouse is on hovered on a point.
    """
    vis = annot.get_visible()
    if event.inaxes == ax:
        cont, ind = scatter.contains(event)
        if cont:
            update_annot(scatter, annot, ind, norm,
                         cmap, labels, mesh_ids, classes)
            annot.set_visible(True)
            fig.canvas.draw_idle()
        else:
            if vis:
                annot.set_visible(False)
                fig.canvas.draw_idle()
                

#Preparation for scatterplot
labels = get_labels()
legend_keys = sorted({int(k):v for k,v in labels.items()}.keys())
legend_labels = []
for key in legend_keys:
    if labels[str(key)] not in legend_labels:
        legend_labels.append(labels[str(key)])
unique_labels = set(sorted(labels.values()))
class_indices = list(range(0, len(unique_labels)))
class_labels = dict(zip(unique_labels, class_indices))

feature = []
indices_classes = []
mesh_ids = []


for index, row in features.iterrows():
    
    mesh_class = labels[str(row.iloc[-1])]
    indices_classes.append(class_labels[mesh_class])
    mesh_ids.append(row.iloc[-1])
    

classes = get_labels()


nr_classes = len(set(classes))
len_classes = len(set(classes.values()))
norm = mpl_colors.Normalize(0, nr_classes)
fig, ax = plt.subplots()
colors = plt.cm.jet(np.linspace(0, 1, len_classes))

#Preparation for dr

features_dr = features.drop(['mesh_id','name'], axis=1)

#dimension reduction
tsne = manifold.TSNE(init="pca", perplexity=21, learning_rate=180)
feature_tsne = tsne.fit_transform(features_dr)

X= feature_tsne
mesh_ids = [int(i) for i in mesh_ids]
#mesh_ids.sort()
scatter = ax.scatter(x=X[:, 0], y=X[:, 1], c= mesh_ids, cmap=matplotlib.colors.ListedColormap(colors),
                         alpha=1, edgecolors='none', s=20)

annot = ax.annotate("", xy=(0, 0), xytext=(10, 10),
                        textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
annot.set_visible(False)

#For legend
recs=[]
for i in range(0,len(set(classes.values()))):
    recs.append(mpatches.Rectangle((0,0),1,1,fc=colors[i]))

ax.legend(recs, legend_labels, bbox_to_anchor=(1, 1.02), loc='upper left')

ax.set_title('Scatterplot of T-SNE procedure')

fig.canvas. mpl_connect("motion_notify_event", lambda event: hover(
        fig, scatter, annot, ax, event, norm, matplotlib.colors.ListedColormap(colors), class_labels, mesh_ids, classes))
plt.show()
