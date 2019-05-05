import numpy as np
import matplotlib.pyplot as plt

# telescope = ['NOT','LT','P60','P200','LCOGT1M','APO','Keck1','WHT']
# colors = ["r", "blue", "#cb416b",'#ff028d','C','y', 'g','M'] 
'''
class Handler(object):
    def __init__(self, color):
        self.color=color
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        patch = plt.Rectangle([x0, y0], width, height, facecolor='blue',
                                   edgecolor='k', transform=handlebox.get_transform())
        patch2 = plt.Rectangle([x0+width/2., y0], width/2., height, facecolor=self.color,
                                   edgecolor='k', transform=handlebox.get_transform())
        handlebox.add_artist(patch)
        handlebox.add_artist(patch2)
        return patch


plt.gca()

telescope = ['NOT','LT','P60','P200','LCOGT1M','APO','Keck1','WHT']
handles = [plt.Rectangle((0,0),1,1) for i  in range(8)]
colors = ["r", "blue", "#cb416b",'#ff028d','C','y', 'g','M'] 
hmap = dict(zip(handles, [Handler(color) for color in colors] ))

for i in telescope:
    plt.legend(handles=handles, labels=str+telescope[i], handler_map=hmap)

plt.show()'''



import matplotlib.pyplot as plt
telescopes = ['NOT','LT','P60','P200','LCOGT1M','ARC','Keck1','WHT','UH88','SOAR']
colors = ["r", "blue", "#e17701",'#ff028d','C','y', 'g','M','#7CFC00','#9932CC'] 
f = lambda m,c: plt.plot([],[],marker=m, color=c, ls="none")[0]
handles = [f("s", colors[i]) for i in range(len(colors))]
labels = telescopes

legend = plt.legend(handles, labels, loc='center', framealpha=1, frameon=False, ncol=2, fontsize='xx-large',markerscale=3)

def export_legend(legend, filename="telescope.png"):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)

export_legend(legend)
#plt.title('Color which representing which telescope')
plt.savefig('telescope.pdf')
plt.show()

