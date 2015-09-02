# -*- coding: utf-8 -*-

from mpl_toolkits.mplot3d import proj3d
import matplotlib.pyplot as plt

from plasma_model_operations_cuda import get_component, p_prev_time, make_dir
from constant import *


def make_3d_plot_with_speed(ax, position, time, n, plot_type='SPEED', annotated=False):
    ax.set_xlim3d([0, n[0]])
    ax.set_ylim3d([0, n[1]])
    ax.set_zlim3d([0, n[2]])
    grid_dim = (int(n[0]/X_STEP), int(n[1]/Y_STEP), int(n[2]/Z_STEP))
    n = position.shape
    labels, flags = None, None
    if annotated:
        labels = [[[False for i in range(grid_dim[0])] for j in range(grid_dim[1])] for k in range(grid_dim[2])]
    flags = np.zeros(grid_dim)
    for num in range(n[1]):
        x_big, y_big, z_big = \
            get_component(position[time][num])
        i, j, k = int(x_big/X_STEP), int(y_big/Y_STEP), int(z_big/Z_STEP)
        ax.scatter(x_big, y_big, z_big, color='red')
        if not annotated:
            continue
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
            if annotated:
                labels[i][j][k] = \
                    ax.annotate(str_label, xy=(x2, y2), xytext=(-15, 15),
                    textcoords='offset points', ha='right', va='bottom',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

    return ax


def make_2d_plot_with_tracks(ax, position, time, n, thinning=1):
    ax.set_xlim([0, n[0]])
    ax.set_ylim([0, n[1]])
    grid_dim = (int(n[0]/X_STEP), int(n[1]/Y_STEP), int(n[2]/Z_STEP))
    n = position.shape
    for num in range(0, n[1], thinning):
        x_big, y_big = \
            get_component(position[time][num], n=2)
        ax.scatter(x_big, y_big, color='red')
        # xes = [position[t][num][0] for t in range(time)]
        # yes = [position[t][num][1] for t in range(time)]
        # ax.plot(xes, yes)

    return ax

def make_tracks_plot_file(prefix, position, time):
    directory = make_dir(prefix, 'by_time')
    n = (X_DIMENSION_GRID, Y_DIMENSION_GRID, Z_DIMENSION_GRID)
    fig = plt.figure()
    ax = fig.add_subplot(211, projection='3d')
    ax2 = fig.add_subplot(212)
    plt.title('time = {} '.format(time))
    ax = make_3d_plot_with_speed(ax, position, p_prev_time(time), n, plot_type='COORD')
    ax2 = make_2d_plot_with_tracks(ax2, position, p_prev_time(time), n)
    # plt.show()
    plt.savefig("{}/carbon_two_coord_{:04d}.png".format(directory, time))
    plt.clf()
    plt.close()

    plt.plot([Mnuc*np.sqrt(position[p_prev_time(time)][i][5]**2 + position[p_prev_time(time)][i][6]**2 + position[p_prev_time(time)][i][7]**2)  for i in range(position.shape[1])])
    plt.title('carbon speed from particle number time = {} '.format(time))
    plt.savefig("{}/carbon_speed_by_pn_time={:04d}".format(directory, time))
    plt.clf()
    plt.close()


def make_tracks_plot_file_with_crashes(prefix, position, crashes, time, thinning=1):
    directory = make_dir(prefix, 'by_time')
    n = (X_DIMENSION_GRID, Y_DIMENSION_GRID, Z_DIMENSION_GRID)
    electron_crashes, helium_crashes, carbon_crashes = crashes
    fig = plt.figure()
    ax = fig.add_subplot(211, projection='3d')
    ax2 = fig.add_subplot(212)
    plt.title('time = {} '.format(time))
    ax = make_3d_plot_with_speed(ax, position, p_prev_time(time), n, plot_type='COORD')
    ax2 = make_2d_plot_with_tracks(ax2, position, p_prev_time(time), n, thinning=thinning)

    xes = [num[0] for num in electron_crashes]
    yes = [num[1] for num in electron_crashes]
    zes = [num[2] for num in electron_crashes]
    if xes:
        ax.scatter(xes, yes, zes, color='green', s=40, marker='x', label='electron')
        ax2.scatter(xes, yes, color='green', s=40, marker='x', label='electron')
    xes = [num[0] for num in helium_crashes]
    yes = [num[1] for num in helium_crashes]
    zes = [num[2] for num in helium_crashes]
    if xes:
        ax.scatter(xes, yes, zes, color='blue', s=35, marker='+', label='helium')
        ax2.scatter(xes, yes, color='blue', s=35, marker='+', label='helium')
    xes = [num[0] for num in carbon_crashes]
    yes = [num[1] for num in carbon_crashes]
    zes = [num[2] for num in carbon_crashes]
    if xes:
        ax.scatter(xes, yes, zes, color='black', s=35, marker='*', label='carbon')
        ax2.scatter(xes, yes, color='black', s=35, marker='*', label='carbon')
    if len(electron_crashes) + len(helium_crashes) + len(carbon_crashes) > 0:
        ax2.legend()
    plt.savefig("{}/carbon_two_coord_{:04d}.png".format(directory, time))
    plt.clf()
    plt.close()

    # plt.plot([Mnuc*np.sqrt(position[p_prev_time(time)][i][5]**2 + position[p_prev_time(time)][i][6]**2 + position[p_prev_time(time)][i][7]**2)  for i in range(position.shape[1])])
    # plt.title('carbon speed from particle number time = {} '.format(time))
    # plt.savefig("{}/carbon_speed_by_pn_time={:04d}".format(directory, time))
    # plt.clf()
    # plt.close()

def make_tracks_plot_file_with_crashes_only_2d(prefix, position, crashes, listen_particles, time, thinning=1):
    directory = make_dir(prefix, 'by_time')
    n = (X_DIMENSION_GRID, Y_DIMENSION_GRID, Z_DIMENSION_GRID)
    electron_crashes, helium_crashes, carbon_crashes = crashes
    fig = plt.figure()
    ax2 = fig.add_subplot(111)
    plt.title('time = {} '.format(time))
    ax2 = make_2d_plot_with_tracks(ax2, position, p_prev_time(time), n, thinning=thinning)
    for num in listen_particles:
        x_big, y_big = \
            get_component(position[p_prev_time(time)][num], n=2)
        ax2.scatter(x_big, y_big, color='green')

    xes = [num[0] for num in electron_crashes]
    # xes = [num[0] for num in carbon_crashes]
    yes = [num[1] for num in electron_crashes]
    # yes = [num[1] for num in carbon_crashes]
    zes = [num[2] for num in electron_crashes]
    # zes = [num[2] for num in carbon_crashes]
    if xes:
        ax2.scatter(xes, yes, color='green', s=40, marker='x', label='electron')
    xes = [num[0] for num in helium_crashes]
    yes = [num[1] for num in helium_crashes]
    zes = [num[2] for num in helium_crashes]
    if xes:
        ax2.scatter(xes, yes, color='blue', s=35, marker='+', label='helium')
    xes = [num[0] for num in carbon_crashes]
    yes = [num[1] for num in carbon_crashes]
    zes = [num[2] for num in carbon_crashes]
    if xes:
        ax2.scatter(xes, yes, color='black', s=35, marker='*', label='carbon')
    if len(electron_crashes) + len(helium_crashes) + len(carbon_crashes) > 0:
        ax2.legend()
    plt.savefig("{}/carbon_two_coord_{:04d}.png".format(directory, time))
    plt.clf()
    plt.close()

    # plt.plot([Mnuc*np.sqrt(position[p_prev_time(time)][i][5]**2 + position[p_prev_time(time)][i][6]**2 + position[p_prev_time(time)][i][7]**2)  for i in range(position.shape[1])])
    # plt.title('carbon speed from particle number time = {} '.format(time))
    # plt.savefig("{}/carbon_speed_by_pn_time={:04d}".format(directory, time))
    # plt.clf()
    # plt.close()


def make_inten_plot_file(prefix, inten, time):
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

def make_distribution_plot_file(prefix, begin, end):
    directory = make_dir(prefix, 'graphs')
    # fig, (ax1, ax2) = plt.subplots(2)
    # ax1.hist(begin, normed=True, histtype='stepfilled')
    # ax2.hist(end, normed=True, histtype='stepfilled')
    # ax1.set_title("distribution at begin")
    # ax2.set_title("distribution at end")
    # fig, ax1 = plt.subplots(1)
    plt.hist(begin, alpha=0.6, normed=True, histtype='stepfilled', label='begin')
    plt.hist(end, alpha=0.6, normed=True, histtype='stepfilled', label='end')
    plt.legend()
    # ax1.set_title("distribution at begin")
    # ax2.set_title("distribution at end")
    plt.savefig("{}/speed_distribution".format(directory))
    plt.clf()
    plt.close()


def make_collision_distribution_plot_file(prefix, typeI, typeII, typeIII):
    directory = make_dir(prefix, 'graphs')
    # fig, (ax1, ax2) = plt.subplots(2)
    # ax1.hist(begin, normed=True, histtype='stepfilled')
    # ax2.hist(end, normed=True, histtype='stepfilled')
    # ax1.set_title("distribution at begin")
    # ax2.set_title("distribution at end")
    # fig, ax1 = plt.subplots(1)
    if typeI != []:
    	plt.hist(typeI, alpha=0.6, normed=True, histtype='stepfilled', label='I')
    if typeII != []:
    	plt.hist(typeII, alpha=0.6, normed=True, histtype='stepfilled', label='II')
    if typeIII != []:
    	plt.hist(typeIII, alpha=0.6, normed=True, histtype='stepfilled', label='III')
    plt.legend()
    # ax1.set_title("distribution at begin")
    # ax2.set_title("distribution at end")
    plt.savefig("{}/collision_distribution".format(directory))
    plt.clf()
    plt.close()


def make_speed_position_plot_file(prefix, plot_data, particles, time):
    directory = make_dir(prefix, 'graphs')
    for pds, p in zip(plot_data, particles):
        fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
        for pd in pds:
            ax1.plot([i[0] for i in pd])
            ax2.plot([i[1] for i in pd])
            ax3.plot([i[2] for i in pd])
        ax1.set_title("x from time particle num={}".format(p))
        ax2.set_title("y from time particle num={}".format(p))
        ax3.set_title("z from time particle num={}".format(p))
        plt.savefig("{}/positions_from_time_pn={}".format(directory, p))
        plt.clf()

        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)
        for pd in pds:
            ax1.plot([i[0] for i in pd], [i[5]*Mnuc if np.abs(i[5]*Mnuc)>1E-30 else 0 for i in pd])#, color=colors[p])
            ax2.plot([i[0] for i in pd], [i[6]*Mnuc if np.abs(i[6]*Mnuc)>1E-30 else 0 for i in pd])#, color=colors[p])
            ax3.plot([i[0] for i in pd], [i[7]*Mnuc if np.abs(i[7]*Mnuc)>1E-30 else 0 for i in pd])#, color=colors[p])
            ax4.plot([i[0] for i in pd], [np.sqrt(i[5]**2 + i[6]**2 + i[7]**2)*Mnuc for i in pd])#, color=colors[p])
        ax1.set_title("v_x from x particle num={}".format(p))
        ax2.set_title("v_y from x particle num={}".format(p))
        ax3.set_title("v_z from x particle num={}".format(p))
        ax4.set_title("|v| from x particle num={}".format(p))
        plt.savefig("{}/speeds_from_x_pn={}.png".format(directory, p))
        plt.clf()

        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)
        for pd in pds:
            ax1.plot([i for i in range(len(pd))], [i[5]*Mnuc if np.abs(i[5]*Mnuc)>1E-30 else 0 for i in pd])#, color=colors[p])
            ax2.plot([i for i in range(len(pd))], [i[6]*Mnuc if np.abs(i[6]*Mnuc)>1E-30 else 0 for i in pd])#, color=colors[p])
            ax3.plot([i for i in range(len(pd))], [i[7]*Mnuc if np.abs(i[7]*Mnuc)>1E-30 else 0 for i in pd])#, color=colors[p])
            ax4.plot([i for i in range(len(pd))], [np.sqrt(i[5]**2 + i[6]**2 + i[7]**2)*Mnuc for i in pd])#, color=colors[p])
        ax1.set_title("v_x from time particle num={}".format(p))
        ax2.set_title("v_y from time particle num={}".format(p))
        ax3.set_title("v_z from time particle num={}".format(p))
        ax4.set_title("|v| from time particle num={}".format(p))
        plt.savefig("{}/speeds_50_50_from_time_pn={}.png".format(directory, p))
        plt.clf()

        fig, (ax1, ax2) = plt.subplots(2, sharex=True)
        for pd in pds:
            ax1.plot([i[0] for i in pd], [i[1] for i in pd])
            ax2.plot([i[0] for i in pd], [i[2] for i in pd])
        ax1.set_title("y from x particle num={}".format(p))
        ax2.set_title("z from x particle num={}".format(p))
        plt.savefig("{}/positions(position)_50_50_pn={}".format(directory, p))

        fig, (ax1) = plt.subplots(1)
        for pd in pds:
            ax1.plot([i[1] for i in pd], [i[2] for i in pd])
        ax1.set_title("z from y particle num={}".format(p))
        plt.savefig("{}/positions(position)_2_50_50_pn={}".format(directory, p))
        plt.clf()
        plt.close()

def make_crash_track(prefix, bp1, ep1, bp2, ep2, r1, r2, vs, time_step, time):
    directory = make_dir(prefix, 'crash')
    t0 = time_step
    x00, y00, z00 = bp1#0, 0, -10
    x01, y01, z01 = ep1#10, 10, 10
    x10, y10, z10 = bp2#0, 0, 10
    x11, y11, z11 = ep2#11, 10, -10
    k0x = (x01 - x00) / t0
    k0y = (y01 - y00) / t0
    k0z = (z01 - z00) / t0
    k1z = (z11 - z10) / t0
    k1y = (y11 - y10) / t0
    k1x = (x11 - x10) / t0
    b0x, b0y, b0z = x00, y00, z00# - k0x*t0
    b1x, b1y, b1z = x10, y10, z10
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter([x00, x01], [y00, y01], [z00, z01], color='blue')
    ax.plot([x00, x01], [y00, y01], [z00, z01], color='blue')
    ax.plot([x00+r1, x00-r1], [y00, y00], [z00, z00], color='blue')
    ax.plot([x00, x00], [y00, y00], [z00+r1, z00-r1], color='blue')
    ax.plot([x00, x00], [y00+r1, y00-r1], [z00, z00], color='blue')
    ax.plot([x01+r1, x01-r1], [y01, y01], [z01, z01], color='blue')
    ax.plot([x01, x01], [y01, y01], [z01+r1, z01-r1], color='blue')
    ax.plot([x01, x01], [y01+r1, y01-r1], [z01, z01], color='blue')
    ax.scatter([x10, x11], [y10, y11], [z10, z11], color='green')
    ax.plot([x10, x11], [y10, y11], [z10, z11], color='green')
    ax.plot([x10-r2, x10+r2], [y10, y10], [z10, z10], color='green')
    ax.plot([x10, x10], [y10-r2, y10+r2], [z10, z10], color='green')
    ax.plot([x10, x10], [y10, y10], [z10-r2, z10+r2], color='green')
    ax.plot([x11-r2, x11+r2], [y11, y11], [z11, z11], color='green')
    ax.plot([x11, x11], [y11-r2, y11+r2], [z11, z11], color='green')
    ax.plot([x11, x11], [y11, y11], [z11-r2, z11+r2], color='green')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    # ax.set_xlim(0, 10)
    # ax.set_ylim(0, 10)
    # ax.set_zlim(0, 10)
    v = vs[0]
    if v != None:
        ax.scatter([k0x * v + b0x, k1x * v + b1x], [k0y * v + b0y, k1y * v + b1y], [k0z * v + b0z, k1z * v + b1z], color='yellow')
    v = vs[1]
    if v != None:
        ax.scatter([k0x * v + b0x, k1x * v + b1x], [k0y * v + b0y, k1y * v + b1y], [k0z * v + b0z, k1z * v + b1z], color='red')
    plt.savefig("{}/crash_{:04d}".format(directory, time))
    plt.clf()
    plt.close()