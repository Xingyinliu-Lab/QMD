# -*- coding: utf-8 -*-
import json


def save_dict(filename, dic):
    '''save dict into json file'''
    with open(filename, 'w') as json_file:
        json.dump(dic, json_file, ensure_ascii=True)


def load_dict(filename):
    '''load dict from json file'''
    with open(filename, "r") as json_file:
        dic = json.load(json_file)
    return dic
