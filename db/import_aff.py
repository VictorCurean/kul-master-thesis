import psycopg2
import csv
from os import listdir
from os.path import isfile, join

DB_NAME = "kul_thesis"
USER = "postgres"
PASSWORD = "admin"
HOST = "localhost"

def get_conn():
    return psycopg2.connect(host=HOST, database=DB_NAME, user=USER, password=PASSWORD)


def check_files():
    path = "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\datasets_aff\\dataset_agg\\"
    files_length = len([f for f in listdir(path) if isfile(join(path, f))])


    summary_file = "C:\\Users\\curea\\Documents\\Thesis Code\\kul-master-thesis\\kul-master-thesis\\kul-thesis-ab\\datasets_aff\\Summary.tsv"
    with open(summary_file, "r") as f:
        lines = len(f.readlines())

    a = 1

