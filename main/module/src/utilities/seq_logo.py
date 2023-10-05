import matplotlib.pyplot as plt
import logomaker as lm
import numpy as np


def get_seq_logo(seq_list, title, save_seq_logo_file):
    df_pfm = lm.alignment_to_matrix(seq_list) # position frequency matrix
    df_ppm = (df_pfm/len(seq_list)) # position probability matrix
    df_pwm = np.log2(df_ppm/0.25) # position weight matrix
    df_ppm = df_ppm*2
    logo = lm.Logo(df_ppm, font_name = 'Arial Rounded MT Bold')
    logo.ax.set_xlabel('Position',fontsize=14)
    logo.ax.set_ylabel('Probability * 2',fontsize=14)
    logo.ax.set_xticks(range(len(seq_list[0])), range(1, len(seq_list[0])+1))
    logo.ax.axhline(1, color='k', linewidth=1, linestyle=':')
    logo.ax.axhline(0.5, color='k', linewidth=1, linestyle=':')
    logo.ax.axhline(1.5, color='k', linewidth=1, linestyle=':')
    plt.title(title)
    plt.savefig(save_seq_logo_file)
    df_pwm['max_column'] = df_pwm.apply(lambda x: x.idxmax(), axis=1)
    canonical_seq = ''.join(df_pwm['max_column'])
    return canonical_seq


