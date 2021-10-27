

import plotly.graph_objects as go
#from peaksIdentification.postprocessing import sliding_avg
import numpy as np


def plotPeaksShifts(peaks, new_peaks, df_results, name1, name2, alg_name='SDS', acc=0.):

    # TUTTI I PICCHI FREE:
    free_xy = peaks[['x', 'y']].to_numpy()
    free_keys = peaks.index.tolist()

    # TUTTI I PICCHI WITH LIGAND:
    with_ligand_xy = new_peaks[['x', 'y']].to_numpy()
    with_ligand_keys = new_peaks.index.tolist()


    #fig = go.Figure()
    fig = go.Figure(
        layout= go.Layout(title="["+alg_name+" (Acc: {:.2f}%)] ".format(acc*100)+name1+" - "+name2,
                          font=dict(size = 10))
    )

    fig.update_layout(width=1300, height=1000)
    fig['layout']['xaxis']['autorange'] = "reversed"
    fig['layout']['yaxis']['autorange'] = "reversed"

    fig.add_trace(
        go.Scatter(
            mode='markers+text',
            x=free_xy[:, 0],
            y=free_xy[:, 1]*0.2,
            marker=dict(
                color='LightSkyBlue',
                size=4,
                line=dict(
                    color='MediumPurple',
                    width=1
                )
            ),
            name='Free Protein Peaks',
            text=free_keys, textposition="bottom center"
        ))

    fig.add_trace(
        go.Scatter(
            mode='markers+text',
            x=with_ligand_xy[:, 0],
            y=with_ligand_xy[:, 1]*0.2,
            marker=dict(
                color='Coral',
                size=4,
                line=dict(
                    color='MediumPurple',
                    width=1
                )
            ),
            name='Protein + Ligand Peaks',
            text=with_ligand_keys, textposition="bottom center"
        ))


    for Pi_key in df_results.index.to_list():
        assigned_to = df_results.loc[Pi_key, 'assigned_to']
        if assigned_to != 'na':
            free_peak = df_results.loc[Pi_key, ['x', 'y']].to_numpy()
            pred_peak = df_results.loc[Pi_key, ['pred_x', 'pred_y']].to_numpy()
            if Pi_key == assigned_to:
                color = "MediumPurple"
            else:
                color = "red"

        fig.add_trace(
            go.Scatter(
                x=[free_peak[0], pred_peak[0]],
                y=[free_peak[1]*0.2, pred_peak[1]*0.2],
                mode='lines',
                showlegend=False,
                text='provaaa',
                line=dict(color=color)))
    fig.show()



def plotProfile(df_results, name1, name2, alg_name='SDS', acc=None):
    estimated_dist_dict = {}
    real_dist_dict = {}

    #print(df_results)
    #print(df_results.columns)
    import sys
    #sys.exit()

    for residue_idx in df_results.index.to_list():
        #print("=======> ", residue_idx)
        estimated_shift = df_results.loc[residue_idx, 'est_shift']
        real_shift = df_results.loc[residue_idx, 'real_shift']


        if estimated_shift > -1.:
            estimated_dist_dict[residue_idx] = estimated_shift
        else:
            estimated_dist_dict[residue_idx] = 0.


        if real_shift > -1.:
            real_dist_dict[residue_idx] = real_shift
        else:
            real_dist_dict[residue_idx] = 0.

    if acc != None:
        fig2 = go.Figure(
            layout= go.Layout(title="["+alg_name+" (Acc: {:.2f}%)] ".format(acc*100)+name1+" - "+name2,
                              font=dict(size = 10)))
    else:
        fig2 = go.Figure(
            layout= go.Layout(title="["+alg_name+"] "+name1+" - "+name2,
                              font=dict(size = 10)))

    fig2.add_trace(go.Bar(
        #x=df1.index,
        #y=df1["dist"],
        x=list(estimated_dist_dict.keys()),
        y=list(estimated_dist_dict.values()),
        name='Estimated Distance',
        marker_color='blue'
    ))

    #fig2 = px.bar(df1, x="Name", y="Real_dist")
    fig2.add_trace(go.Bar(
        #x=df1.index,
        x = list(real_dist_dict.keys()),
        y=list(real_dist_dict.values()),
        name='Real Distance',
        marker_color='lightgray'
    ))

    fig2.show()


def plotHistogram(df1, new_peaks, name1="name1", name2="name2", alg_name="algorithm", accuracy=-1, divide_one_fifth = False):

    real_dist_dict = None

    if divide_one_fifth:
        df1['y'] = df1['y'] / 5.
        new_peaks['y'] = new_peaks['y'] / 5.

    estimated_dist_dict = {}

    for residue_idx in df1.index:
        assigned_res_idx = df1.loc[residue_idx, 'assigned_to']
        if assigned_res_idx != 'na':
            jjj = df1.loc[residue_idx, ['x', 'y']].to_numpy()
            iii = new_peaks.loc[assigned_res_idx, ['x', 'y']].to_numpy()
            kk = np.sqrt(((iii - jjj)**2).sum())
            estimated_dist_dict[residue_idx] = kk
        else:
            estimated_dist_dict[residue_idx] = 0.

    real_dist_dict = {}
    for residue_idx in df1.index:
        if residue_idx in new_peaks.index:
            jjj = df1.loc[residue_idx, ['x', 'y']].to_numpy()
            iii = new_peaks.loc[residue_idx, ['x', 'y']].to_numpy()
            kk = np.sqrt(((iii - jjj)**2).sum())
            real_dist_dict[residue_idx] = kk

    #print(real_dist_dict)
    #print(list(real_dist_dict.keys()))

    #print("REAL DIST DICT ", real_dist_dict) # real distance dictionary
    #print("SIZE REAL DIST DICT ", len(real_dist_dict))
    #print("= = = = =>\n",df1)

    #df1[ df1['assigned_to']=='na' ]['dist'] = 0.
    # df1.loc[df1['assigned_to']=='na', 'dist'] = 0.
    distances = df1['dist'].tolist()

    #sl_avg = sliding_avg(distances, half_window_size=3) # sliding window avg on our estimated shift distances
    #print(sl_avg)

    #print(len(df1['Name'].tolist()), df1['Name'].tolist())
    #print(len(df1['Index'].tolist()), df1['Index'].tolist())
    #print("==============")


    #df1['window_avg'] = sl_avg
    #sl_avg = df1['window_avg'].tolist()

    # plot dello istogramma
    #fig2 = px.bar(df1, x = "Name", y="Distance")
    fig2 = go.Figure(
        layout= go.Layout(title="["+alg_name+"] "+name1+" - "+name2+"({:.2f})".format(accuracy*100),
                          font=dict(size = 10))
    )


    fig2.add_trace(go.Bar(
        #x=df1.index,
        #y=df1["dist"],
        x=list(estimated_dist_dict.keys()),
        y=list(estimated_dist_dict.values()),
        name='Estimated Distance',
        marker_color='blue'
    ))


    #fig2 = px.bar(df1, x="Name", y="Real_dist")
    fig2.add_trace(go.Bar(
        #x=df1.index,
        x = list(real_dist_dict.keys()),
        y=list(real_dist_dict.values()),
        name='Real Distance',
        marker_color='lightgray'
    ))

    # plotta il punto indicante il valore della media mobile
    '''fig2.add_trace(
        go.Scatter(
            mode='markers+text+lines',
            x=df1.index,
            y=sl_avg,
            marker=dict(
                color='Coral',
                size=4,
                line=dict(color='MediumPurple',width=1
                )
            ),
            name='Window avg', text='', textposition="bottom center"
        ))'''

    fig2.show()
    #return sl_avg, df1['Name'].tolist(), df1['Index'].tolist()