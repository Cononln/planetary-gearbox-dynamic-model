"""Fault-dependent kinematic row lengths for signed vibration processing."""
import math
import numpy as np


def recurrence(fault, sun=21, planet=31, ring=84, planets=3):
    teeth = sun if fault == 'sun' else planet if fault == 'planet' else None
    if teeth is None:
        raise ValueError('Only sun and planet recurrence are defined')
    turns = teeth//math.gcd(teeth, ring)
    mesh = turns*ring
    events = mesh*planets//sun if fault == 'sun' else mesh//planet
    return dict(carrier_turns=turns, mesh_cycles=mesh, fault_cycles=events)


def angular_rows(motion, repeat, points, ring=84, start=3., stop=59.):
    carrier = motion['phase']/(2*np.pi*ring)
    assert np.all(np.diff(carrier) > 0)
    origin = np.interp(start, motion['time'], carrier)
    n = int(np.floor((np.interp(stop, motion['time'], carrier)-origin)/repeat))
    if n < 4:
        raise ValueError('Insufficient complete rows for the sun pilot')
    edges_angle = origin+np.arange(n+1)*repeat
    edges = np.interp(edges_angle, carrier, motion['time'])
    within = np.arange(points)/points*repeat
    query = np.interp((edges_angle[:-1, None]+within).ravel(), carrier, motion['time']).reshape(n, points)
    return edges, query


def split_row_correlation(rows):
    x = np.asarray(rows, float)
    even = x[::2].mean(axis=0)
    odd = x[1::2].mean(axis=0)
    even -= even.mean()
    odd -= odd.mean()
    den = np.linalg.norm(even)*np.linalg.norm(odd)
    return float(np.dot(even, odd)/den) if den > 1e-20 else float('nan')
