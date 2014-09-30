import numpy as np
import h5py
import MDAnalysis

def tail(f, n, offset=0):
  command = "tail -n {} {}".format(n+offset, f)
  stdin,stdout = os.popen2(command)
  stdin.close()
  lines = stdout.readlines(); stdout.close()
  return lines[0:-offset]

def copygro(gro_in, gro_out, frames):
    with open(gro_in) as f:
        i = 0
        for l in f:
            l_arr = l.split()
            if "trjconv" in l_arr:
                t0 = float(l_arr[-1])
                print "t0 = {}".format(t0)
                atoms = int(f.next())
                print "atoms = {}".format(atoms)
                break
        for atom in xrange(atoms):
            f.next()
        box = f.next().strip()
        l1 = f.next()
        t1 = float(l1.split()[-1])
        print "box = {}".format(box)
        print "t1 = {}".format(t1)
        print "dt = {}".format(t1-t0)

    lf = tail(gro_in, 1, atoms+2)[0]
    tf = float(lf.split()[-1])
    print "tf = {}".format(tf)
    N  = int(round((tf / (t1-t0)))) + 1

    if np.any(np.array(frames) >= N):
        raise ValueError("Frame value too high in copygro, frames must be less than {}".format(N))
    
    f_out = open(gro_out,'w')
    for frame in frames:
        f_out = open(gro_out,'a')
        with open(gro_in) as f:
            n = 0
            while n <= frame:
                print 'Frame ', n, "/", N
                l0 = f.next()
                l1 = f.next()
                if n == frame:
                    f_out.write(l0)
                    f_out.write(l1)
                for atom in xrange(atoms):
                    l = f.next()
                    if n == frame:
                        f_out.write(l)
                lf = f.next()
                if n == frame:
                    f_out.write(lf)
                n += 1
        f_out.close()


def trr2hdf5(trr_in, hdf_out, hdf_atom):
    print "Loading universe..."
    universe=MDAnalysis.Universe(trr_in[1], trr_in[0])
    print "Getting the dt step..."
    print universe.trajectory.dt

def gro2hdf5(gro_in, hdf_out, hdf_atom):
    with open(gro_in) as f:
        i = 0
        for l in f:
            l_arr = l.split()
            if "trjconv" in l_arr:
                t0 = float(l_arr[-1])
                print "t0 = {}".format(t0)
                atoms = int(f.next())
                print "atoms = {}".format(atoms)
                break
        for atom in xrange(atoms):
            f.next()
        box = f.next().strip()
        l1 = f.next()
        t1 = float(l1.split()[-1])
        print "box = {}".format(box)
        print "t1 = {}".format(t1)
        print "dt = {}".format(t1-t0)

    lf = tail(gro_in, 1, atoms+2)[0]
    tf = float(lf.split()[-1])
    print "tf = {}".format(tf)
    N  = int(round((tf / (t1-t0)))) + 1

    h5 = h5py.File(hdf_out,'w')
    ds = h5.create_dataset(hdf_atom, (N,atoms,3), dtype='f32')
    with open(gro_in) as f:
        for n in xrange(N):
            print 'Frame ', n, "/", N
            l0 = f.next()
            l1 = f.next()
            for atom in xrange(atoms):
                if atom % 100000 == 0:
                    print "\tatom",atom
                l = f.next()
                pos = l[20:].strip()
                pos = np.array([float(xyz) for xyz in pos.split()])
                ds[n,atom,:] = pos[:]
            lf = f.next()
    h5.close()
