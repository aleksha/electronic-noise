#===============================================================================
# CODE MUST BE REFACTORED !!!
#===============================================================================
# see for details
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.fftpack.rfft.html#scipy.fftpack.rfft
#===============================================================================
dump_file_path = "./dump.txt"
N_EVENTS       = 100
N_CHANNELS     = 2692
BAN_LEVEL      = 65
#===============================================================================
import shelve
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab   as mlab
from   scipy.fftpack import  fft,  ifft
from   scipy.fftpack import rfft, irfft
from   scipy.stats   import norm, chi2, beta
from math import pow, sqrt, atan, pi
#===============================================================================
def Arg( x, y ):
    if x > 0             : return atan(y/x)
    if x < 0  and y >= 0 : return atan(y/x)+pi
    if x < 0  and y < 0  : return atan(y/x)-pi
    if x == 0 and y > 0  : return 0.5*pi
    if x == 0 and y == 0 : return 0            # not defined actually
#===============================================================================
class anode_noise:
    """Class, which containg info about noise using calibration sample
       obtained with pulse generator"""

    def __init__(self, dump_path):
        """Load csv dump, make fft"""
        self.channels = N_CHANNELS
        self.events   = 0
        print("Loading anode dump file")
        self.path = dump_path
        self.dataset = []

        self.model = []
        self.is_model = False

        max_diff  = []
        down_diff = []
        ban_list  = []
        ev_num    = 0
        dumpfile = open(self.path,"r")
        for line in dumpfile:
            value_str_list = (line[:-1]).split(",")
            value_list = []
            for idx in range(self.channels):
                value_list.append( float( value_str_list[idx] ) )
            max_level = max( sum(value_list)/len(value_list) - min(value_list) ,
                             max(value_list) - sum(value_list)/len(value_list) )
            max_diff .append ( max_level )
            if  sum(value_list)/len(value_list) - min(value_list) > max(value_list) - sum(value_list)/len(value_list)   :
                down_diff.append ( sum(value_list)/len(value_list) - min(value_list)  )
            if max_level > BAN_LEVEL:
                ban_list.append( ev_num )
            ev_num = ev_num + 1
        dumpfile.close()
        n, bins, patches = plt.hist( np.array(max_diff) , 100 , facecolor='green', alpha = 0.5)
        m, bbbb, patche2 = plt.hist( np.array(down_diff), bins = bins, facecolor='yellow', alpha = 0.5)
        plt.grid( True )
        plt.xlabel("maximal difference from event base line")
        plt.ylabel("events")
        plt.savefig( "MAX_DIFF.png" )
        plt.clf()

        dumpfile = open(self.path,"r")
        ev_num = 0
        for line in dumpfile:
            if ev_num not in ban_list :
                value_str_list = (line[:-1]).split(",")
                value_list = []
                for idx in range(self.channels):
                    value_list.append( float( value_str_list[idx] ) )
                self.dataset .append( np.array( value_list ) )
                self.events   = self.events + 1
            ev_num = ev_num + 1
        dumpfile.close()

        print("CHANNELS : " + str( self.channels ) )
        print("EVENTS   : " + str( self.events   ) )
        fft_list  = []
        fft_list2 = []
        print("Performing fast Fourier transformation")
        for event in self.dataset:
            fft_list .append(  fft( event ) )
        self.fftset  = fft_list

    def draw_event(self, num, fig_name = "EVENT.png"):
        """Draw event number num into figure fig_name"""
        plt.plot( self.dataset[num] )
        plt.grid( True )
        plt.xlabel("time, ch.")
        plt.ylabel("val, a.u.")
        plt.title(r'Signal for EVENT %d' %(num))
        plt.savefig( fig_name )
        print("Time spectrum for event " + str(num) + " into " + fig_name )
        plt.clf()
        return self.dataset[num]

    def draw_extended_event(self, num, fig_name = "EXTD_EVENT.png"):
        """Draw event number num into figure fig_name"""
        evt_fft = np.copy ( self.fftset[num] )
        evt_newl = []

        for ii in range(0,self.channels):
            evt_newl.append( evt_fft[ii]*2. )
            evt_newl.append( 0 )

        evt_new = np.array(evt_newl)
        ext_event  = ifft( evt_new )

        plt.plot( ext_event, "r" )
        plt.plot( self.dataset[num] )
        plt.grid( True )
        plt.xlabel("time, ch.")
        plt.ylabel("val, a.u.")
        plt.title(r'Signal for EVENT %d' %(num))
        plt.savefig( fig_name )
        print("Time spectrum for event " + str(num) + " into " + fig_name )
        plt.clf()

    def filter_event(self, num, flt = "gap 100 250", draw_sig = True, fig_name = "FLTR_EVENT.png", xmin=0, xmax=2692):
        """Draw event number num into figure fig_name"""
        all_flt = ["gap","cut"]
        filter = flt.split(" ")
        if filter[0] in all_flt:
            evt_fft = np.copy ( self.fftset[num] )
            evt_newl_part1 = []
            evt_newl_part2 = []
            if filter[0]=="gap":
                evt_newl_part1.append( evt_fft[0] )
                for ii in range(1,1+len(evt_fft)/2):
                    if ii>int(filter[1]) and ii<int(filter[2]):
                        evt_newl_part1.append(   evt_fft[ii]             )
                        evt_newl_part2.insert(0, evt_fft[ii].conjugate() )
                    else:
                        evt_newl_part1.append(   0. )
                        evt_newl_part2.insert(0, 0. )
            elif filter[0]=="cut":
                evt_newl_part1.append( evt_fft[0] )
                for ii in range(1,1+len(evt_fft)/2):
                    if ii>int(filter[1]) and ii<int(filter[2]):
                        evt_newl_part1.append(   0. )
                        evt_newl_part2.insert(0, 0. )
                    else:
                        evt_newl_part1.append(   evt_fft[ii]             )
                        evt_newl_part2.insert(0, evt_fft[ii].conjugate() )

            evt_newl = evt_newl_part1 + evt_newl_part2
#            for ii in range(len(evt_fft)):
#                print(str(ii)+"\t"+ str(evt_fft[ii]) + "\t" + str(evt_newl[ii]) )

            evt_new = np.array(evt_newl)
            ext_event  = ifft( evt_new )

            plt.plot( ext_event, "r" )
            axes = plt.gca()
            axes.set_xlim([xmin,xmax])
            if draw_sig:
                plt.plot( self.dataset[num] )
            plt.grid( True )
            plt.xlabel("time, ch.")
            plt.ylabel("val, a.u.")
            plt.title("Signal for EVENT " + str(num) + " with FILTER " + flt)
            plt.savefig( fig_name )
            print("Time spectrum for event " + str(num) + " into " + fig_name )
            plt.clf()

    def filter_average(self, ev_min , ev_max, flt = "none", fig_name = "FLTR_AVE.png",
                             xmin=0, xmax=2692, fmin=0, fmax = 200, freq_name = "FLTR_SPEC.png"):
        """Draw event number num into figure fig_name"""
        condition = ev_min<ev_max and ev_max < self.events and ev_min>-1
        if condition :
            all_flt = ["gap","cut","none"]
            filter = flt.split(" ")
            if filter[0] in all_flt:
                base = 0.
                nevt = 0.
                for ch in range( len( self.dataset[0] ) ):
                    for ev in range( ev_min, ev_max ):
                        base = base + float( self.dataset[ev][ch] )
                        nevt = nevt + 1.
                base = base / nevt

                signal = []
                for ch in range( len( self.dataset[0] ) ):
                    val = 0.
                    for ev in range( ev_min, ev_max ):
                        val = val + float( self.dataset[ev][ch] )
                    signal.append( val / float( ev_max - ev_min ) - base )
#                if ch < 300:
#                    print(str(ch) + "\t" + str( val / float( len( self.dataset ) ) ))

                sig = np.array( signal )
                plt.plot( sig, "green" )
                axes = plt.gca()
                axes.set_xlim([xmin,xmax])
                plt.grid( True )
                plt.xlabel("time, ch.")
                plt.ylabel("val, a.u.")
                plt.title( "Signal for events FROM " + str(ev_min) + " TO " + str(ev_max) +" with FILTER " + flt )
                plt.savefig( fig_name )
                print("Time spectrum for events into " + fig_name )
                plt.clf()

                spec = fft( sig ) / float( self.channels )
                plt.plot( np.abs(spec[fmin:fmax]), "green" )
                plt.grid( True )
                plt.xlabel("freq.")
                plt.ylabel("val, a.u.")
                plt.title( "Spectrum for events FROM " + str(ev_min) + " TO " + str(ev_max) +" with FILTER " + flt )
                plt.savefig( freq_name )
                print("Time spectrum for events into " + freq_name )
                plt.clf()

                return sig

    def draw_spectrum(self, num, fig_name = "SPECTRUM.png"):
        """Draw absolute value of frequiency spectrum"""
        plt.plot( np.abs( self.fftset[num][1:] ) )
        plt.grid( True )
        plt.xlabel("freq., 1/ch.")
        plt.ylabel("val, a.u.")
        plt.title(r'Spectrum for EVENT %d' %(num))
        plt.savefig( fig_name )
        print("Distribution of absolute value ( ev = " + str(num) +  " ) to " + fig_name )
        print(" (zero term substructed)")
        print(" len ( freq ) = " + str( len( self.fftset[0] ) ) )
        plt.clf()


    def create_model( self ):
        """Re / Im gaussian fits for all channels"""

        for chan in range( int( 1 + self.channels/2 ) ):
            #print(" MODEL FOR CHANNEL :" + str(chan))

            data_re  = []
            data_im  = []

            for evt_fft in self.fftset:
                data_re .append( np.real ( evt_fft[chan] )/self.channels )
                data_im .append( np.imag ( evt_fft[chan] )/self.channels )

            mu1, si1  = norm.fit(data_re)
            if chan != 0 or chan != self.channels/2:
                mu2, si2  = norm.fit(data_im)
            else:
                mu2 = 0
                si2 = 0

            self.model.append( (mu1, si1, mu2, si2) )

        self.is_model = True

    def generate_event( self , fig_name = "GEN_EVENT.png"):
        """Return random event"""
        if not self.model:
            self.create_model()

        gen_list = []

        evt_newl_part1 = []
        evt_newl_part2 = []
        for chan in range( int( 1 + self.channels/2 ) ):

            mu1 = self.model[chan][0]
            si1 = self.model[chan][1]
            mu2 = self.model[chan][2]
            si2 = self.model[chan][3]

            rep = norm(mu1,si1).rvs()
            if mu2==0 and si2==0:
                imp = 0.0
            else :
                imp = norm(mu2,si2).rvs()
            evt_newl_part1.append(    rep + 1j*imp )
            if chan!=0:
                evt_newl_part2.insert( 0, rep - 1j*imp )

        evt_newl = evt_newl_part1 + evt_newl_part2
        evt_new = np.array(evt_newl) * self.channels
        gen_event  = ifft( evt_new )

#        plt.plot( gen_event, "black" )
#        plt.grid( True )
#        plt.xlabel("time, ch.")
#        plt.ylabel("val, a.u.")
#        plt.title("Generated EVENT")
#        plt.savefig( fig_name )
#        print("Time spectrum for generated event stored into " + fig_name )
#        plt.clf()
        return gen_event

    def create_model_chan(self, chan, fig_name = "CHANNEL.png", save_fig = True, test_it = False):
        """Re / Im gaussian fits for chanel"""
        data_abs = []
        data_arg = []
        data_re  = []
        data_im  = []
        for evt_fft in self.fftset:
            data_abs.append( np.abs  ( evt_fft[chan] )/self.channels )
            data_arg.append( np.angle( evt_fft[chan] )               )
            data_re .append( np.real ( evt_fft[chan] )/self.channels )
            data_im .append( np.imag ( evt_fft[chan] )/self.channels )


        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2) #, sharex='col', sharey='row')
        f.suptitle('FREQ ' + str(chan) )
        n1, b1, p1 = ax1.hist(data_re , 50, facecolor='green', normed = True ,label = "real part"      )
        n2, b2, p2 = ax2.hist(data_im , 50, facecolor='green', normed = True ,label = "imaginary part" )
        n3, b3, p3 = ax3.hist(data_abs, 50, facecolor='green', normed = True ,label = "absolute value" )
        n4, b4, p4 = ax4.hist(data_arg, 50, facecolor='green', normed = True ,label = "complex phase"  )

        ax1.legend(loc='best')
        ax2.legend(loc='best')
        ax3.legend(loc='best')
        ax4.legend(loc='best')

        mu1, si1  = norm.fit(data_re)
        lnspc1    = np.linspace(min(data_re), max(data_re), len(data_re))
        pdf_norm1 = norm.pdf(lnspc1, mu1, si1)
        ax1.plot(lnspc1, pdf_norm1, "r", linewidth = 3)

        if test_it:
            re_arr = norm(mu1,si1).rvs(size = 20*len(data_abs))

        mu2, si2  = norm.fit(data_im)
        lnspc2    = np.linspace(min(data_im), max(data_im), len(data_im))
        pdf_norm2 = norm.pdf(lnspc2, mu2, si2)
        ax2.plot(lnspc2, pdf_norm2, "r", linewidth = 3)
        if test_it:
            im_arr = norm(mu2,si2).rvs(size = 20*len(data_abs))
            abs_arr = np.sqrt(np.power(re_arr,2) + np.power(im_arr,2))
            arg_arr = np.angle(re_arr + 1j*im_arr)

            ax3.hist(abs_arr, bins = b3, facecolor="yellow", normed = True , alpha = 0.5)
            ax4.hist(arg_arr, bins = b4, facecolor="yellow", normed = True , alpha = 0.5)

        if save_fig :
            plt.savefig( fig_name )
        plt.clf()



    def draw_average_spectrum(self, fmin = 1, fmax = 100,  fig_name = "SPECTRUM_AVE.png"):
        """Draw absolute value of frequiency spectrum"""

        means = []
        for chan in range(fmin,fmax):
            datos = []
            for evt_fft in self.fftset :
                datos.append( np.abs( evt_fft[chan] ) )
            (mu, sigma) = norm.fit(datos)
            means.append( mu/self.channels )
        plt.plot( means )
        plt.grid( True )
        plt.xlabel("freq., 1/num.of.ch.")
        plt.ylabel("mean abs. val., a.u.")
        plt.title('Average spectrum')
        plt.savefig( fig_name )
        print("Distribution of absolute value to " + fig_name )
        print("  in range = ("+str(fmin) + ","+str(fmax) + ")")
        plt.clf()

    def hist_abs(self, chan, fig_name = "ABS_VAL.png"):
        """Distribution of abs. val. for frequency channel for all events"""
        datos = []
        for evt_fft in self.fftset:
            datos.append( np.abs( evt_fft[chan] )/self.channels )
        n, bins, patches = plt.hist(datos, 100, facecolor='green', normed = True)

        xt = plt.xticks()[0]
        xmin, xmax = min(xt), max(xt)
        lnspc = np.linspace(xmin, xmax, len(datos))
        ab,bb,cb,db = beta.fit(datos)
        pdf_beta = beta.pdf(lnspc, ab, bb, cb, db)
        plt.plot(lnspc, pdf_beta,  "r", label="Beta", linewidth = 3)

        plt.grid( True )
        plt.xlabel("abs value")
        plt.ylabel("events")
        plt.title(r'Distribution for FREQ %d' %(chan))
        plt.savefig( fig_name )
        print("Distribution of absolute values for all events ( freq. = " + str(chan) + " ) saved into " + fig_name )
        plt.clf()

    def hist_arg(self, chan, fig_name = "ARG_VAL.png"):
        """Distribution of phase for frequency channel for all events"""
        datos = []
        for evt_fft in self.fftset:
            datos.append( np.angle( evt_fft[chan] ) )
        n, bins, patches = plt.hist(datos, 50, facecolor='green', normed = True)

#        xt = plt.xticks()[0]
#        xmin, xmax = min(xt), max(xt)
#        lnspc = np.linspace(xmin, xmax, len(datos))
#        hist = np.histogram(dataos, bins=50)
#        hist_dist = rv_histogram(hist)
#        plt.plot(lnspc, hist_dist.pdf(lnspc), "r", label="Beta", linewidth = 3)

        plt.grid( True )
        plt.xlabel("complex phase")
        plt.ylabel("events")
        plt.title(r'Distribution for FREQ %d' %(chan))
        plt.savefig( fig_name )
        print("Distribution of compex phase for all events ( freq. = " + str(chan) + " ) saved into " + fig_name )
        plt.clf()

    def check_args(self, chi2_alpha = 0.997, fig_name = "CHCK_ARGS.png"):
        """Make chi2 test that argument distributions are flat"""
        print("\n\n====== CHECK THAT PHASE IS FLAT ======\n")
        for chan in range(1,self.channels/2):
            datos = []
            for evt_fft in self.fftset:
                datos.append( np.angle( evt_fft[chan] ) )
            n, bins, patches = plt.hist(datos, self.events/100, facecolor='green')
            chi2_obs = 0.
            for n_i in n:
                chi2_obs = chi2_obs + float( pow( n_i-100 ,2) )/float(n_i)
            if chi2_obs > chi2.interval(chi2_alpha,len(n))[1]:
                print("-----> CHANNEL " + str(chan) + "\t OBS. :" + str(chi2_obs) + "\t EXP.: " + str(chi2.interval(chi2_alpha,len(n)) ) )
            elif chi2_obs < chi2.interval(chi2_alpha,len(n))[0]:
                print("-----> CHANNEL " + str(chan) + "\t OBS. :" + str(chi2_obs) + "\t TOO GOOD TO BE TRUE" )
            else:
                print("-----> CHANNEL " + str(chan) + "\t OBS. :" + str(chi2_obs) + "\t FINE" )
            plt.clf()
        print("\n======      END PHASE CHECK     ======\n")

    def check_args_corr(self, chan_min   = 0 , chan_stop=256, chan_max = 256,
                              chan_depth = 10, chan2_shift = 0,
                              chi2_alpha = 0.997, fig_name = "CHCK_DIFF.png"):
        """Make chi2 test that chan1 chan2 argument difference distributions are flat"""
        print("\n\n====== CHECK THAT DIFF. IS FLAT ======\n")

        rdb = shelve.open("temp_store.shelve")
        tmp_arr_2d = rdb["phase_diff_array"]

        if len(tmp_arr_2d)<chan_max:
            arr_2d = np.zeros( (chan_max,chan_max) )
            for ch1 in range(0,chan_max):
                for ch2 in range(0,chan_max):
                    arr_2d[ch1][ch2] = 0.5
            for ch1 in range(len(tmp_arr_2d)):
                for ch2 in range(0,len(tmp_arr_2d)):
                    arr_2d[ch1][ch2] = tmp_arr_2d[ch1][ch2]
        else:
            arr_2d = tmp_arr_2d

        for chan1 in range(chan_min,chan_stop):
            print("-----> CHANNEL " + str(chan1) + " / " + str(chan_stop) )
            for chan2 in range(chan1+chan2_shift,chan1+chan2_shift+chan_depth):
                if chan2 < chan_max and chan1<chan_max:
                    if arr_2d[chan1][chan2]>0. or arr_2d[chan1][chan2]<1.:
                        datos = []
                        for evt_fft in self.fftset:
                            if np.angle( evt_fft[chan1] ) - np.angle( evt_fft[chan2] ) > 0 :
                                datos.append( np.angle( evt_fft[chan1] ) - np.angle( evt_fft[chan2] ) )
                            else :
                                datos.append( np.angle( evt_fft[chan1] ) - np.angle( evt_fft[chan2] ) + pi*2. )
                        n, bins, patches = plt.hist(datos, self.events/100, facecolor='green')
                        plt.title(r'Difference of phases between %d and %d' %(chan1,chan2))
                        plt.savefig("TEMP.png")
                        if chan1 == chan2 :
                            arr_2d[chan1][chan2]=0.
                        else:
                            chi2_obs = 0.
                            for n_i in n:
                                chi2_obs = chi2_obs + float( pow( n_i-100 ,2) )/float(n_i)
                            if chi2_obs > chi2.interval(chi2_alpha,len(n))[1]:
                                arr_2d[chan1][chan2]=1.
                            else:
                                arr_2d[chan1][chan2]=0.
                        plt.clf()

        for ch1 in range(0,chan_max):
            for ch2 in range(ch1+1,chan_max):
                arr_2d[ch2][ch1] = arr_2d[ch1][ch2]

        plt.imshow(arr_2d,cmap="gray", interpolation='none')
        plt.tight_layout()
        plt.xlabel("channel")
        plt.ylabel("channel")
        plt.title(r'Non-flat phases difference between channels')
        plt.savefig(fig_name)
        plt.clf()
        print("\n======      END DIFF. CHECK     ======\n")
        rdb.close()
        return arr_2d

    def corr_scatter(self, chan1, chan2, fig_name = "CORR_SCT.png"):
        """Scatter plot for two channels"""
        datos1 = []
        datos2 = []
        datos3 = []
        datos4 = []
        for evt_fft in self.fftset:
            datos1.append( np.real( evt_fft[chan1] ) )
            datos2.append( np.imag( evt_fft[chan1] ) )
            datos3.append( np.real( evt_fft[chan2] ) )
            datos4.append( np.imag( evt_fft[chan2] ) )
        plt.scatter(datos1, datos2, s=1, c = "red" , alpha=0.5)
        plt.scatter(datos3, datos4, s=1, c = "blue", alpha=0.5)
        plt.grid( True )
        plt.savefig( fig_name )
        plt.clf()

    def corr_abs(self, chan1, chan2, fig_name = "ABS_CORR.png"):
        """Scatter plot for two complex phases"""
        datos1 = []
        datos2 = []
        for evt_fft in self.fftset:
            datos1.append( np.abs( evt_fft[chan1] ) )
            datos2.append( np.abs( evt_fft[chan2] ) )
        plt.scatter(datos1, datos2, s=10, c = "red", alpha=0.5)
        plt.grid( True )
        plt.xlabel("abs. val. for freq. " + str(chan1))
        plt.ylabel("abs. val. for freq. " + str(chan2))
        plt.title(r'Correlation of abs. val. between %d and %d' %(chan1,chan2))
        plt.savefig( fig_name )
        print("Corrwlation for all events ( freq. = " + str(chan1) + " and " + str(chan2) + " ) saved into " + fig_name )
        plt.clf()

    def phase_abs(self, chan1,fig_name = "ABS_ARG.png"):
        """Scatter plot for abs and complex phases"""
        datos1 = []
        datos2 = []
        for evt_fft in self.fftset:
            datos1.append( np.angle( evt_fft[chan1] ) )
            datos2.append( np.abs  ( evt_fft[chan1] ) )
        plt.scatter(datos1, datos2, s=10, c = "green", alpha=0.5)
        plt.grid( True )
        plt.savefig( fig_name )
        plt.clf()

    def corr_arg(self, chan1, chan2, fig_name = "ARG_CORR.png"):
        """Scatter plot for two complex phases"""
        datos1 = []
        datos2 = []
        for evt_fft in self.fftset:
            datos1.append( np.angle( evt_fft[chan1] ) )
            datos2.append( np.angle( evt_fft[chan2] ) )
        plt.scatter(datos1, datos2, s=10, c = "red", alpha=0.5)
        plt.grid( True )
        plt.xlabel("phase for freq. " + str(chan1))
        plt.ylabel("phase for freq. " + str(chan2))
        plt.title(r'Correlation of phases between %d and %d' %(chan1,chan2))
        plt.savefig( fig_name )
        print("Corrwlation for all events ( freq. = " + str(chan1) + " and " + str(chan2) + " ) saved into " + fig_name )
        plt.clf()

    def diff_arg(self, chan1, chan2, fig_name = "ARG_DIFF.png"):
        """Diggerence for two complex phases"""
        datos = []
        for evt_fft in self.fftset:
            if np.angle( evt_fft[chan1] ) - np.angle( evt_fft[chan2] ) > 0 :
                datos.append( np.angle( evt_fft[chan1] ) - np.angle( evt_fft[chan2] ) )
            else :
                datos.append( np.angle( evt_fft[chan1] ) - np.angle( evt_fft[chan2] ) + pi*2. )
        n, bins, patches = plt.hist(datos, 30, facecolor='red', alpha = 0.75)
        plt.grid( True )
        plt.xlabel("phase " + str(chan1) + " - phase " + str(chan2))
        plt.ylabel("events")
        plt.title(r'Difference of phases between %d and %d' %(chan1,chan2))
        plt.savefig( fig_name )
        print("Difference for all events ( freq. = " + str(chan1) + " and " + str(chan2) + " ) saved into " + fig_name )
        plt.clf()

#===============================================================================
anode = anode_noise(dump_file_path)

evt = anode.draw_event(1024)
sig = anode.filter_average( ev_min =    0, ev_max = 4850, flt="none" ,fmax = 200)

plt.plot( evt - ( sum(evt) / len(evt) ) , "blue"  )
plt.plot( sig, "green" )
axes = plt.gca()
plt.grid( True )
plt.xlabel("time, ch.")
plt.ylabel("val, a.u.")
plt.title( "Signal for one event and averaged one")
plt.savefig( "EV_EV-AVE.png" )
plt.clf()

#anode.filter_average( ev_min = 2000, ev_max = 4000, flt="none" )
#anode.filter_average( flt="none" , xmin =0, xmax = 300 )
#anode.filter_event(1024,flt="cut 28 30",xmin =0, xmax = 1000)
#anode.filter_event(1020,flt="gap 100 250", draw_sig = False)
#anode.draw_extended_event(20)
#anode.draw_spectrum(20)
#anode.draw_average_spectrum(1,200)
#anode.hist_abs(3)
#anode.hist_arg(3)
#anode.phase_abs(29)
#anode.diff_arg(100,150)
#anode.corr_arg(100,150)
#anode.corr_abs(100,150)
#anode.corr_scatter(100,100)

#anode.check_args()

#a2d = anode.check_args_corr(chan_min = 1, chan_stop = 2, chan_max = 400, chan_depth = 1, chan2_shift = 0)
#rdb = shelve.open("temp_store.shelve")
#rdb["phase_diff_array"] = a2d
#rdb.close()

#for idx in range(1,8):
#    if idx<10:
#        sname = "0"+str(idx)
#    else:
#        sname = str(idx)
#    anode.hist_arg(idx,fig_name = "ARG_" + sname + ".png")

#anode.create_model_chan( 110 , test_it = True)

from statistics import mean, stdev

out_file = open("noise_events.data","w")
for gev in range(N_EVENTS):
    real_arr = np.real ( anode.generate_event() )
#    print(len(real_arr))
    mean_arr = mean(real_arr)
    print(gev,mean_arr,stdev(real_arr),sep="\t")
    ss = ""
    for tch in real_arr:
        ss += "{:4.2f}".format(tch-mean_arr) + " "
    ss += "\n"
    out_file.write(ss)
out_file.close()

