#!/usr/bin/env python3
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("f1", help="First input file for SE/PE.", type=str)
    parser.add_argument("-f2", help="Second input file for PE.", type=str, default=None)
    parser.add_argument("-a1", help="First adapter for SE/PE adapter trimming", type=str, default=None)
    parser.add_argument("-a2", help="Second adapter for SE/PE adapter trimming", type=str, default=None)
    parser.add_argument("-e", help="maximum error rate", type=float, default=None)
    parser.add_argument("-t", help="thread num", type=int, default=1)
    parser.add_argument("-m", help="min length for SE/PE adapter trimming", type=int, default=0)
    parser.add_argument("-idx", help="the indexing table of peat", type=str)
    parser.add_argument("-mp", help="Earrings single end seed length", type=int, default=50)
    parser.add_argument("-mr", help="Earrings PE reverse mismatch rate", type=float, default=0.7)
    parser.add_argument("-sc", help="Earrings PE sequence compare rate", type=float, default=0.9)
    parser.add_argument("-ac", help="Earrings PE adapter compare rate", type=float, default=0.8)
    parser.add_argument("-pf", help="Earrings prune factor", type=float, default=0.03)
    parser.add_argument("-se", help="Earrings sensitive mode", type=int, default=0)
    parser.add_argument("-af", help="atropos additional flag", type=str)
    parser.add_argument("-da", help="disable given adapters", type=int, default=1)
    parser.add_argument("-px", help="enable polyX trimming in fastp", type=int, default=0)
    parser.add_argument("-rc", help="Record memory usage", type=int, default=0)

    return parser.parse_args()
