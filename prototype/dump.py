# -*- coding: utf-8 -*-
import time
import numpy as np
import os.path
import pickle
from constant import CONSTANT_VALUES
from functools import reduce
import sqlite3 as lite
from datetime import date as dt

__author__ = 'fedul0x'


def is_dump(filename):
    try:
        with open(filename, 'rb') as dump_file:
            _ = pickle.load(dump_file)
    except:
        return False
    return True


def restore_from_dump(dumpfolder, conds):
    if dumpfolder is None:
        return None
    files = map(lambda x: os.path.join(dumpfolder, x), os.listdir(dumpfolder))
    files = filter(lambda x: is_dump(x), files)
    for file in files:
        dump_file = open(file, 'rb')
        consts = pickle.load(dump_file)

        def check(x):
            if conds.get(x[0], None) is None or\
            x[0] == 'TIME_STEP':
                return True
            else:
                return x[1] == conds[x[0]]
        # check = lambda x: x[1] == conds[x[0]]
        # print(reduce(lambda x, y: x and y, (map(check, consts.items()))))
        if reduce(lambda x, y: x and y, (map(check, consts.items()))):
            return pickle.load(dump_file)
    return None


def save_to_dump(dumpfolder, conds, data):
    if dumpfolder is None:
        return None
    filename = 'dump_{}.pickle'.format(time.strftime('%Y%m%d%H%M%S'))
    filename = os.path.join(dumpfolder, filename)
    with open(filename, 'wb') as dump_file:
        pickle.dump(conds, dump_file)
        pickle.dump(data, dump_file)


class DbConnection:

    def __init__(self, dbfile, descrpt, date=None):
        self.connection = lite.connect(dbfile)
        self.cursor = self.connection.cursor()
        date = dt.today() if date is None else date
        sql = '''insert into experiment (date, description)
            VALUES ("{}", "{}")'''.format(date, descrpt)
        self.cursor.execute(sql)
        self.connection.commit()
        self.experiment_last_id = \
            self.cursor.execute('select max(id) from experiment').fetchone()[0]
        # self.collision_last_id = \
        #     self.cursor.execute('select max(id) from collision').fetchone()[0]
        self.particle_last_id = \
            self.cursor.execute('select max(id) from particle').fetchone()[0]
        self.iteration_last_id = \
            self.cursor.execute('select max(id) from iteration').fetchone()[0]

    def new_iteration(self, time):
        sql = '''insert into iteration (id_experiment, time)
            VALUES ("{}", "{}")'''.format(self.experiment_last_id, time)
        self.cursor.execute(sql)
        self.connection.commit()
        self.iteration_last_id = self.cursor.execute('select max(id) from iteration').fetchone()[0]

    def new_particle(self, carbon, tension, currtime, carbonum):
        for i in range(carbonum):
            c = ', '.join([str(carbon[currtime][i][j]) for j in range(9)])
            t = ', '.join([str(tension[i][j]) for j in range(3)])
            sql = '''insert into particle (id_iteration, pos_x, pos_y, pos_z,
                radius, charge, speed_x, speed_y, speed_z, guid,
                tension_x, tension_y, tension_z)
                VALUES ({}, {}, {})'''\
                .format(self.iteration_last_id, c, t)
            self.cursor.execute(sql)
        self.connection.commit()
        sql = 'select max(id) from particle'
        self.particle_last_id = self.cursor.execute(sql).fetchone()[0]

    def new_collision(self, collisions):
        for i in collisions:
            # TODO rewrite i[0]+self.particle_last_id, i[1]+self.particle_last_id
            # TODO calc energy
            id1, id2, e = \
                i[0]+self.particle_last_id, i[1]+self.particle_last_id, i[2]
            sql = '''insert into collision (id_iteration,
                id_particle_1, id_particle_2, energy)
                VALUES ({}, {}, {}, {})'''\
                .format(self.iteration_last_id, id1, id2, e)
            self.cursor.execute(sql)
        self.connection.commit()

    def new_final(self, carbons, tension):
        for i in range(len(carbons)):
            c = ', '.join([str(carbons[i][j]) for j in range(9)])
            t = ', '.join([str(tension[i][j]) for j in range(3)])
            sql = '''insert into final (id_iteration, pos_x, pos_y, pos_z,
                radius, charge, speed_x, speed_y, speed_z, guid,
                tension_x, tension_y, tension_z)
                VALUES ({}, {}, {})'''\
                .format(self.iteration_last_id, c, t)
            self.cursor.execute(sql)
        self.connection.commit()
