#%% Import modules
import baltic as bt
import matplotlib as mpl
from matplotlib import pyplot as plt
%matplotlib inline

#%% Parse input (trees & attributes)
trees = {  'int': bt.loadNexus('/Users/Sidney/Dropbox/siv/beast/discreteTraits/h24/asymmetric/original/trees/mcc/6_raw.mcc', absoluteTime=False),
            'env': bt.loadNexus('/Users/Sidney/Dropbox/siv/beast/discreteTraits/h24/asymmetric/original/trees/mcc/9_raw.mcc', absoluteTime=False)}
# 'gag': bt.loadNexus('/Users/Sidney/Dropbox/siv/beast/discreteTraits/h24/asymmetric/original/trees/mcc/2_raw.mcc', absoluteTime=False),

def parse_attribs(path, asfloat=False):
    splitlines = [line.strip().split('\t') for line in open(path, 'r') if len(line.split('\t')) == 2]
    if asfloat:
        return {k:float(v) for k,v in splitlines}
    else:
        return {k:v.strip() for k,v in splitlines}

# common_names = parse_attribs('/Users/Sidney/Dropbox/siv/figures/common_names.tsv')
colors = parse_attribs('/Users/Sidney/Dropbox/siv/figures/colors.tsv')


#%% Collect genome coordinates & breakpoints
#     'pol': [2085, 5096, 1]
CDS = {'vpu': [6045, 6062, 2], 'tat': [8379, 8469, 3], "3' LTR": [9086, 9719, 2], 'vif': [5041, 5619, 3], 'vpr': [5559, 5850, 1], 'int': [4438, 5096, 1], 'rev': [8379, 8653, 3], "5' LTR": [1, 634, 3], 'nef': [8797, 9168, 3], 'pol': [2085, 4438, 1], 'env': [6225, 8795, 1], 'gag': [790, 2292, 3]}

#%%

breakpoints = [1343.,1951.,2291.,3526.,4438.,5558.,5896.,6567.,6974.,7716.,8616.]

three_prime = ['tat', 'vpu', 'vpr', 'rev', 'env', 'rev', 'nef', "3' LTR"]
middle = ['vif', 'int']


#%% Plot genome map & breakpoints
fig,ax = plt.subplots(figsize=(50,10),facecolor='w')
genomeL=9719.0

def rescale(x):
    return (float(x)/genomeL)

for gene, i in CDS.items():
    length = rescale(i[1]-i[0])
    if gene == 'vpu':
        ax.text(rescale(i[0])+length, 0.1*i[2],'%s'%(gene),va='center',ha='left',size=24,zorder=11)
    elif gene == 'rev':
        ax.text(rescale(i[0])+length, 0.1*i[2],'%s'%(gene),va='center',ha='left',size=24,zorder=11)
    elif gene == 'tat' and i[0] == 5831:
        ax.text(rescale(i[0])-length/1.5, 0.1*i[2],'%s'%(gene),va='center',ha='left',size=24,zorder=11)
    elif gene == 'tat':
        ax.text(rescale(i[0])-length*1.3, 0.1*i[2],'%s'%(gene),va='center',ha='left',size=24,zorder=11)
    elif gene == 'int':
        pass
    else:
        ax.text(rescale(i[0])+0.5*length, 0.1*i[2],'%s'%(gene),va='center',ha='center',size=24,zorder=11)

    if gene in three_prime:
        c = colors['Mona_Monkey']
    elif gene in middle:
        c = colors['Red-capped_Mangabey']
    else:
        c = 'lightgray'

    if gene != 'pol':
        ax.arrow(rescale(i[0]), 0.1*i[2], length, 0.0, alpha=0.6,head_width=0.07, width=0.07,head_length=0.05*rescale(i[1]-i[0]),length_includes_head=True,facecolor=c)
    else:
        ax.arrow(rescale(i[0]), 0.1*i[2], length, 0.0, alpha=0.6,head_width=0.07, width=0.07,head_length=0.0,length_includes_head=True,facecolor=c)

plt.eventplot([rescale(bp) for bp in breakpoints], orientation='horizontal', lineoffsets=0,
          linelengths=1, linewidths=5, linestyles='dashed', color='k')
ax.text(0, -0.05, '0', size=24, ha='center')
ax.text(1, -0.05, '9719', size=24, ha='center')
for bp in breakpoints:
    ax.text(rescale(bp), -0.05, '%d'%int(bp), size=24, ha='center')

#%% Draw trees for segments 2, 6, and 9 in the correct position

cumulative_displace=0 ## this tracks the "current" x position, so trees are plotted one after another

branchWidth=5 ## increase branch width, since trees will be smaller

tree_names=range(ntrees) ## define order in which dict will be accessed

tip_positions={x:{} for x in tree_names} ## remember the position of each tip in each tree

traitName = 'absent'

for tr in tree_names: ## iterate over trees
    cur_tree=trees[tr] ## fetch tree object
    for k in cur_tree.Objects: ## iterate over branches
        if isinstance(k,bt.leaf): ## only interested in leaves
#             print k.name
            tip_positions[tr][k.name]=(k.height,k.y) ## remember tree, tip's position

frac_pos = 0.0

cmap = mpl.cm.viridis
for t,tr in enumerate(tree_names): ## iterate over trees
    cur_tree=trees[tr] ## fetch tree object

    for k in cur_tree.Objects: ## iterate over branches
    #     x=k.x ## or from x position determined earlier
        x=k.height*1.1 ## or use absolute time instead
        y=k.y ## get y position from .drawTree that was run earlier, but could be anything else

    #     xp=k.parent.x ## get x position of current object's parent
        xp=k.parent.height*1.1 ## get x position of current object's parent
        if x==None: ## matplotlib won't plot Nones, like root
            x=0.0
        if xp==None:
            xp=x

        x+=cumulative_displace ## adjust branch position by displacement, which depends on the position of tree in the overall plot
        xp+=cumulative_displace ## same for branch's parent
    #     c='indianred' ## colour can be fixed
    #     c=cmap(k.height/ll.treeHeight) ## or be a function of something else
    #     c=[cmap(k.traits['posterior']) if k.traits.has_key('posterior') else cmap(1.0)][0]
        if k.traits.has_key(traitName):
            c=['indianred' if k.traits[traitName]=='V' else 'steelblue'][0] ## can be discrete too
        else:
            c='k'
        if isinstance(k,bt.leaf): ## if leaf...
            #x=decimalDate(k.name.split('_')[-1],variable=True) ## get x position from name

            s=60 ## tip size can be fixed

            ax.scatter(x,y,s=s,facecolor=cmap(frac_pos),edgecolor='none',zorder=11) ## plot circle for every tip
            ax.scatter(x,y,s=s+0.8*s,facecolor='k',edgecolor='none',zorder=10) ## plot black circle underneath

            if t+1<len(tree_names) and k.name in tip_positions[0] and k.name in tip_positions[tree_names[t+1]]:
                pos_in_first_tree=tip_positions[0][k.name][1] ## fetch y coordinate of same tip in the first tree
                frac_pos=pos_in_first_tree/float(len(cur_tree.Objects))*2.0 ## normalize coordinate to be within interval [0.0,1.0]


                if t!=len(tree_names)-1: ## as long as we're not at the last tree - connect tips with coloured lines
                    next_x,next_y=tip_positions[tree_names[t+1]][k.name] ## fetch coordinates of same tip in next tree
                    next_x+=cumulative_displace+cur_tree.treeHeight+1.5 ## adjust x coordinate by current displacement and future displacement

                    ax.plot([x,next_x],[y,next_y],lw=2,ls='-',color=cmap(frac_pos),zorder=0) ## connect current tip with same tip in the next tree
            else:
                pass

        elif isinstance(k,bt.node): ## if node...
            ax.plot([x,x],[k.children[-1].y,k.children[0].y],lw=branchWidth,color=c,ls='-',zorder=9) ## plot vertical bar

        ax.plot([xp,x],[y,y],lw=branchWidth,color=c,ls='-',zorder=9) ## always plot branch

    cumulative_displace+=cur_tree.treeHeight+1.5 ## increment displacement by the height of the tree




#%% Set figure attributes, save.
ax.set_ylim(0,1) ## set y limits
ax.set_xlim(0, 1)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ax.tick_params(axis='x',labelsize=0,size=0)
ax.tick_params(axis='y',labelsize=0,size=0)
plt.show()

# plt.savefig('sivcpz_2003.png', bbox_inches='tight')
plt.show()
