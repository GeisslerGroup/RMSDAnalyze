import numpy as np

def AtomRMSD(data_tik, output_file, rmsd_step, time_init=0.0, time_step=0.0):
    f = open(output_file,'w') 
    for t0 in xrange(data_tik.shape[0]-rmsd_step):
        f.write("{}\n".format(data_tik.shape[1]))
        f.write("RMSD generated in XYZ by RMSDAnalyze.py @jhaberstroh: t = {}\n".format(time_init + time_step * t0))
        print "outputting time {} of {}".format(t0, data_tik.shape[0]-rmsd_step-1)
        r_ik = data_tik[t0+rmsd_step,:,:] - data_tik[t0,:,:]
        rmsd_i = np.sqrt(np.sum(np.square(r_ik), axis=1))
        rmsd_i[rmsd_i > 5.0] = 0
        #plt.hist(rmsd_i, bins=100)
        #plt.show()
        for rmsd in rmsd_i:
            f.write("C\t{}\t0.0\t0.0\n".format(rmsd))
    f.close()


def AtomActivity(data_tik, output_file, openflag, frames, activity_step, r0, k):
    f = open(output_file, openflag)
    if frames+activity_step > data_tik.shape[0]:
        raise AttributeError("Frames need for activity calculation ({}) exceeds available frames ({})".format(frames+activity_step, data_tik.shape[0]))
    print "Frames: {}".format(frames)
    for t0 in xrange(frames):
        f.write("{}\n".format(data_tik.shape[1]))
        f.write("Activity generated in XYZ by RMSDAnalyze.py @jhaberstroh: t = {}\n".format(frames))
        r_ik = data_tik[t0+activity_step,:,:] - data_tik[t0,:,:]
        rmsd_i = np.sqrt(np.sum(np.square(r_ik), axis=1))

        activity_i = np.reciprocal(1 + np.exp( -k * (rmsd_i - r0) ))
        for activity in activity_i: 
            f.write("C\t{}\t0.0\t0.0\n".format(activity))
    f.close()

