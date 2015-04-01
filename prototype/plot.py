# -*- coding: utf-8 -*-

from mpl_toolkits.mplot3d import proj3d

from plasma_model_operations_cuda import get_component, p_prev_time 
from constant import *


def make_3d_plot_with_speed(ax, position, time, n, plot_type='SPEED'):
    ax.set_xlim3d([0, n[0]])
    ax.set_ylim3d([0, n[1]])
    ax.set_zlim3d([0, n[2]])
    grid_dim = (int(n[0]/X_STEP), int(n[1]/Y_STEP), int(n[2]/Z_STEP))
    n = position.shape
    labels = [[[False for i in range(grid_dim[0])] for j in range(grid_dim[1])] for k in range(grid_dim[2])]
    flags = np.zeros(grid_dim)
    for num in range(n[1]):
        x_big, y_big, z_big = \
            get_component(position[time][num])
        i, j, k = int(x_big/X_STEP), int(y_big/Y_STEP), int(z_big/Z_STEP)
        ax.scatter(x_big, y_big, z_big, color='red')
        if (j == k and not flags[0][j][k]):
            flags[0][j][k] = True
            x2, y2, _ = proj3d.proj_transform(x_big, y_big, z_big, ax.get_proj())
            str_label = ''
            if plot_type == 'SPEED':
                values = get_component(position[time][num], b=5)
                str_label = '{:5.7f},{:5.7f},{:5.7f}'.\
                    format(values[0]*Mnuc, values[1]*Mnuc, values[2]*Mnuc)
            if plot_type == 'COORD':
                values = get_component(position[time][num])
                str_label = '{:5.7f},{:5.7f},{:5.7f}'.\
                    format(values[0], values[1], values[2])
            labels[i][j][k] = \
                ax.annotate(str_label, xy=(x2, y2), xytext=(-15, 15),
                textcoords='offset points', ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

    return ax


def make_2d_plot_with_tracks(ax, position, time, n):
    ax.set_xlim([0, n[0]])
    ax.set_ylim([0, n[1]])
    grid_dim = (int(n[0]/X_STEP), int(n[1]/Y_STEP), int(n[2]/Z_STEP))
    n = position.shape
    labels = [[[False for i in range(grid_dim[0])] for j in range(grid_dim[1])] for k in range(grid_dim[2])]
    for num in range(n[1]):
        x_big, y_big = \
            get_component(position[time][num], n=2)
        ax.scatter(x_big, y_big, color='red')
        xes = [position[t][num][0] for t in range(time)]
        yes = [position[t][num][1] for t in range(time)]
        ax.plot(xes, yes)

    return ax

def make_tracks_plot_file(prefix, position, time):
    directory = make_dir(prefix, 'by_time')
    n = (X_DIMENSION_GRID, Y_DIMENSION_GRID, Z_DIMENSION_GRID)
    fig = plt.figure()
    ax = fig.add_subplot(211, projection='3d')
    ax2 = fig.add_subplot(212)
    plt.title('time = {} '.format(time))
    ax = make_3d_plot_with_speed(ax, positionCarbon, p_prev_time(time), n, plot_type='COORD')
    ax2 = make_2d_plot_with_tracks(ax2, positionCarbon, p_prev_time(time), n)
    # plt.show()
    plt.savefig("{}/carbon_two_coord_{:04d}.png".format(directory, time))
    plt.clf()
    plt.close()

    plt.plot([Mnuc*np.sqrt(positionCarbon[p_prev_time(time)][i][5]**2 + positionCarbon[p_prev_time(time)][i][6]**2 + positionCarbon[p_prev_time(time)][i][7]**2)  for i in range(positionCarbon.shape[1])])
    plt.title('carbon speed from particle number time = {} '.format(time))
    plt.savefig("{}/carbon_speed_by_pn_time={:04d}".format(directory, time))
    plt.clf()
    plt.close()


def make_intent_plot_file(prefix, intent, time):
    directory = make_dir(prefix, 'inten')
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    for i, ax in zip(inten, [ax1, ax2, ax3]):
        ax.plot([k for k,_,_ in i], color='red')
    plt.savefig("{}/int_time={:04d}".format(directory, time))
    plt.clf()
    plt.close()


def make_potential_plot_file(prefix, next_phi, time):
    directory = make_dir(prefix, 'phi')
    n = next_phi.shape
    plot_phi = np.zeros((n[0], n[1]))
    for k in range(n[2]):
        for i in range(n[0]):
            for j in range(n[1]):
                plot_phi[i][j] = next_phi[i][j][k]
        plt.contourf(plot_phi.T)
        plt.colorbar()
        plt.savefig("{}/phi_time={:04d}_z={:02d}".format(directory, time, k))
        plt.clf()
    plt.close()


def make_intensity_plot_file(prefix, intensity, time):
    directory = make_dir(prefix, 'intensity')
    n = intensity.shape;
    plot_intensity = np.zeros((n[0], n[1]))
    for k in range(n[2]):
        for i in range(n[0]):
            for j in range(n[1]):
                # plot_intensity[i][j] = np.sqrt(np.sum([intensity[i][j][k][l]**2 for l in range(3)]))
                plot_intensity[i][j] = intensity[i][j][k][0]

        plt.contourf(plot_intensity.T, cmap=plt.cm.flag)
        plt.colorbar()
        plt.savefig("{}/intensity_time={:04d}_z={:02d}".format(directory, time, k))
        plt.clf()
    plt.close()


def make_tension_plot_file(prefix, tension, time):
    directory = make_dir(prefix, 'tension')
    tension = carbon_tension
    n = tension.shape
    plot_tension = np.zeros(n[0])
    plot_tension_x = np.zeros(n[0])
    plot_tension_y = np.zeros(n[0])
    plot_tension_z = np.zeros(n[0])
    for i in range(n[0]):
        plot_tension[i] = \
            np.sqrt(np.sum([tension[i][l]**2 for l in range(3)]))
        plot_tension_x[i] = tension[i][0]
        plot_tension_y[i] = tension[i][1]
        plot_tension_z[i] = tension[i][2]

    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    ax1.plot(plot_tension_x)
    ax2.plot(plot_tension_y)
    ax3.plot(plot_tension_z)
    plt.savefig("{}/tension_by_x_y_z_time={:04d}".format(directory, time))
    plt.clf()

    plt.plot(plot_tension)
    plt.savefig("{}/tension_by_sum(x,y,z)_time={:04d}".format(directory, time))
    plt.clf()
    plt.close()

def make_distribution_plot_file(prefix, begin, end, time):
    directory = make_dir(prefix, 'graphs')
    fig, (ax1, ax2) = plt.subplots(2)
    ax1.hist(begin_speed_distribution_data, normed=True, histtype='stepfilled')
    ax2.hist(end_speed_distribution_data, normed=True, histtype='stepfilled')
    ax1.set_title("distribution at begin")
    ax2.set_title("distribution at end")
    plt.savefig("{}/speed_distribution".format(directory, p))
    plt.clf()
    plt.close()


def make_speed_position_plot_file(prefix, plot_data, particles, time):
    directory = make_dir(prefix, 'graphs')
    # colors = ['red', 'green', 'blue']
    for pd, p in zip(plot_data, particles):
        fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
        ax1.plot( [i[0] for i in pd])  #, color=colors[p])
        ax1.set_title("x from time particle num={}".format(p))
        ax2.plot( [i[1] for i in pd])#, color=colors[p])
        ax2.set_title("y from time particle num={}".format(p))
        ax3.plot( [i[2] for i in pd])#, color=colors[p])
        ax3.set_title("z from time particle num={}".format(p))
        plt.savefig("{}/positions_50_50_pn={}".format(directory, p))
        plt.clf()

        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)
        ax1.plot([i[0] for i in pd], [i[5]*Mnuc if np.abs(i[5]*Mnuc)>1E-30 else 0 for i in pd])#, color=colors[p])
        ax1.set_title("v_x from x particle num={}".format(p))
        ax2.plot([i[0] for i in pd], [i[6]*Mnuc if np.abs(i[6]*Mnuc)>1E-30 else 0 for i in pd])#, color=colors[p])
        ax2.set_title("v_y from x particle num={}".format(p))
        ax3.plot([i[0] for i in pd], [i[7]*Mnuc if np.abs(i[7]*Mnuc)>1E-30 else 0 for i in pd])#, color=colors[p])
        ax3.set_title("v_z from x particle num={}".format(p))
        ax4.plot([i[0] for i in pd], [np.sqrt(i[5]**2 + i[6]**2 + i[7]**2)*Mnuc for i in pd])#, color=colors[p])
        ax4.set_title("|v| from x particle num={}".format(p))
        plt.savefig("{}/speeds_50_50_pn={}.png".format(directory, p))
        plt.clf()

        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)
        ax1.plot([i for i in range(len(pd))], [i[5]*Mnuc if np.abs(i[5]*Mnuc)>1E-30 else 0 for i in pd])#, color=colors[p])
        ax1.set_title("v_x from time particle num={}".format(p))
        ax2.plot([i for i in range(len(pd))], [i[6]*Mnuc if np.abs(i[6]*Mnuc)>1E-30 else 0 for i in pd])#, color=colors[p])
        ax2.set_title("v_y from time particle num={}".format(p))
        ax3.plot([i for i in range(len(pd))], [i[7]*Mnuc if np.abs(i[7]*Mnuc)>1E-30 else 0 for i in pd])#, color=colors[p])
        ax3.set_title("v_z from time particle num={}".format(p))
        ax4.plot([i for i in range(len(pd))], [np.sqrt(i[5]**2 + i[6]**2 + i[7]**2)*Mnuc for i in pd])#, color=colors[p])
        ax4.set_title("|v| from time particle num={}".format(p))
        plt.savefig("{}/speeds_50_50_from_time_pn={}.png".format(directory, p))
        plt.clf()

        fig, (ax1, ax2) = plt.subplots(2, sharex=True)
        ax1.plot([i[0] for i in pd], [i[1] for i in pd])
        ax1.set_title("y from x particle num={}".format(p))
        ax2.plot([i[0] for i in pd], [i[2] for i in pd])#, color=colors[p])
        ax2.set_title("z from x particle num={}".format(p))
        plt.savefig("{}/positions(position)_50_50_pn={}".format(directory, p))

        fig, (ax1) = plt.subplots(1)
        ax1.plot([i[1] for i in pd], [i[2] for i in pd])  #, color=colors[p])
        ax1.set_title("z from y particle num={}".format(p))
        plt.savefig("{}/positions(position)_2_50_50_pn={}".format(directory, p))
        plt.clf()
        plt.close()