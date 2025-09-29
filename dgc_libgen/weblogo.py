import logomaker
from .utils.seq_mutations import NAT_AAs

def draw_weblogo(weight_mat_df, fix_positions, ax):
    ww_logo = logomaker.Logo(weight_mat_df[NAT_AAs], ax=ax,
                             color_scheme='NajafabadiEtAl2017', vpad=.1, width=.8)

    # style using Logo methods
    ww_logo.style_spines(visible=False)
    ww_logo.style_spines(spines=['left', 'bottom'], visible=True)
    ww_logo.style_xticks(rotation=45, fmt='%d', anchor=0) #

    # style using Axes methods
    ww_logo.ax.set_ylabel("bits", labelpad=-1)
    ww_logo.ax.set_xticklabels(list(weight_mat_df.pos))
    ww_logo.ax.xaxis.set_ticks_position('none')
    ww_logo.ax.xaxis.set_tick_params(pad=-1)
    
    # highlight positions
    for pos in fix_positions:
        ww_logo.highlight_position(p=list(weight_mat_df.pos).index(pos), color='gold', alpha=.5)