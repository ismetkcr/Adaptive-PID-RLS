import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
from create_table_fpdf2 import PDF

class Model:
    def __init__(self, model_type, *args):
        if model_type == 'linear':
            if len(args) != 3:
                raise ValueError("For Linear model, provide three coefficients. (Kp, Tau, Delay)")               
            self.TF = TF(args)
        elif model_type == 'parametric':
            if len(set(len(arg) for arg in args)) !=1:
                raise ValueError("Parameters and model equation matrices should be the same length!")
            self.parametricmodel = ParametricModel(args)
        else:
            raise ValueError("Invalid model type. Choose either 'linear' or 'parametric'.")
    
    

class TF:
    def __init__(self, args):
        self.Kp = args[0]
        self.Tau = args[1]
        self.delay = args[2]
        self.uk_delay = np.zeros(self.delay)
        
    def evaluate(self, dt, yprev, uprev):
        a = -1/self.Tau
        b = self.Kp/self.Tau
        
        ynew = yprev + dt * (a * yprev + b * self.uk_delay[-1])
        self.uk_delay = np.roll(self.uk_delay, 1)
        self.uk_delay[0] = uprev
        return ynew
    
    def dispparameters(self):
        print('Kp = %.f, Tau = %.f, delay = %.f' % (self.Kp, self.Tau, self.delay))
        
    
class ParametricModel:
    def __init__(self, args):
        self.parameters = args[0]
        self.param_list = args[1]
        output_terms = []
        param_dict = {'y1': 'y(k-1)', 'y2': 'y(k-2)', 'u1': 'u(k-1)', 'u2': 'u(k-2)'}
        for param in self.param_list:
          key, power = param.split(',')  # Parametre adını ve üssü ayır
          if key in param_dict:
              if power == '1':
                  output_terms.append(param_dict[key])
              else:
                  output_terms.append(f"{param_dict[key]}^{power}")
        # if output_terms:
        #     print(' + '.join(output_terms)) 
        
        self.model_eqn = output_terms
        
    def cons_model_matrice(self, yprev, yprev2, uprev, uprev2):
        model_matrices = []
        for term in self.model_eqn:
            key, power = term.split('^') if '^' in term else (term, '1')
            power = int(power)
            
            if key == 'y(k-1)':
                model_matrices.append(yprev ** power)
            elif key == 'y(k-2)':
                model_matrices.append(yprev2 ** power)
            elif key == 'u(k-1)':
                model_matrices.append(uprev ** power)
            elif key == 'u(k-2)':
                model_matrices.append(uprev2 ** power)
                
        return np.array(model_matrices)
        
    def evaluate(self, dt, yprev, yprev2, uprev, uprev2):
        model_matrices = self.cons_model_matrice(yprev, yprev2, uprev, uprev2)
        model_output = self.parameters.dot(model_matrices)
        return model_output



class PIDcontroller():
    def __init__(self, dt):
        self.dt = dt
        self.eprev2 = 0
        self.veloutput = 0
    
    def assign_PID_parameters(self, Kc, Ti, Td):
        self.Kc = Kc
        self.Ti = Ti
        self.Td = Td
        
        
    
    def calculate_output(self, r, y, eprev, umin, umax, uprev, s, ubackprev, PID_backprev,  antiwindup = None, tt = None ):
        error = r - y
        P_out = self.Kc * error
        D_out = self.Kc * self.Td * (error - eprev) / self.dt
        
        
        
        if antiwindup is None:
            s += (self.Kc / self.Ti) * error * self.dt
            I_out = s    
            PID_out = P_out + I_out + D_out
            uback = 0
            PID_back = 0
            
            PID_out = max(umin, min(PID_out, umax))
            
            
    
        if antiwindup == 'conditionalintegral':
            if umin < uprev < umax:
                s += (self.Kc / self.Ti) * error * self.dt
            I_out = s    
            PID_out = P_out + I_out + D_out
            uback = 0
            PID_back = 0
            # Clamp output within [umin, umax]
            PID_out = max(umin, min(PID_out, umax))
                
                
    
        if antiwindup == 'backcalculation':
            s += ((1 / tt) * (ubackprev - PID_backprev) + (self.Kc / self.Ti) * error) * self.dt
            I_out = s    
            PID_back = P_out + I_out + D_out
            
            if PID_back <= umin:
                uback = umin
            elif PID_back >= umax:
                uback = umax
            else:
                uback = PID_back 
            
            PID_out = PID_back
            
            PID_out = max(umin, min(PID_out, umax))
            
        return PID_out, error, s, uback, PID_back
    
    def calculate_veloutput(self, r, y, eprev, uprev, umin, umax):
        error = r - y
        a0 = self.Kc * (1 + (self.dt / self.Ti) + (self.Td / self.dt))
        a1 = self.Kc * (-1  - (2 * (self.Td / self.dt)))
        a2 = (self.Kc * self.Td) / self.dt
        self.veloutput += (a0 * error) + (a1 * eprev) + (a2 * self.eprev2)
        
        
        PID_out = self.veloutput
        
        PID_out = max(umin, min(PID_out, umax))
        return PID_out, error, 0, 0, 0
        
    
    def calculate_performance(self, error, iseprev, iaeprev):
        ise = iseprev +  error**2*self.dt
        iae = iaeprev +  np.abs(error)*self.dt
               
        return ise, iae
    
    
        

    
        
class calculatePIDparameters:
    def __init__(self, Kp, Tau, delay, KU, PU):
        self.Kp = Kp
        self.Tau = Tau
        self.delay = delay
        self.KU = KU
        self.PU = PU
        
    def zieglernicholas(self):
        KC = self.KU / 1.7
        TI = self.PU / 2.0
        TD = self.PU / 8.0
    
        return KC, TI, TD
        
    
    def tyreusluyben(self):
        KC = self.KU / 2.2
        TI = self.PU * 2.0
        TD = self.PU / 6.3
        
        return KC, TI, TD
    
    def cohencoon(self):
        KC = (self.Tau/(self.Kp*self.delay)) * ((16 + (3*self.delay/self.Tau))/12)
        TI = (self.delay * (32 + (6*self.delay/self.Tau))) / (13 + (8*self.delay/self.Tau))
        TD = (4*self.delay) / (11 + (2*self.delay/self.Tau))
        
        return KC, TI, TD
        
    
    def itaedist(self):
        KC = (1.357 * (self.delay / self.Tau) ** (-0.947)) / self.Kp
        TI = (self.Tau) / (0.842 * (self.delay / self.Tau) ** (-0.738))
        TD = self.Tau * 0.381 * (self.delay / self.Tau) ** 0.995
        
        return KC, TI, TD


def merge_close_values(arr, threshold):
    merged_values = []
    #current_group = [arr[0]]
    current_group = [0]

    for i in range(1, len(arr)):
        if abs(arr[i] - current_group[-1]) <= threshold:
            current_group.append(arr[i])
        else:
            merged_values.append(np.mean(current_group))
            current_group = [arr[i]]

    # 
    merged_values.append(np.mean(current_group))

    return merged_values


def checkgrad(yf, yfprev, yk, dt, number, t, mtrec=None, ttrec=None, last_peak_value = None, last_peak_time = None, founded_peak = None):
    
    #bu fonksiyona hem filtrelenmiş y değerleri, hem ham y değerleri beslenecek
    #gradyan sıfıra yakınlık ölçümü filtreli veriler üzerinden yapılacak.
    #gradyanın sıfıra yakın olduğu noktalardaki matris dolumu gerçek y değerleri ile hesaplanacak
    #yani  grad = hesabına filtreli y değerlerini socakağım.
    if mtrec is None: #mtrec matrisi, gradyan sıfıra yakın iken gerçek y değerlerini topluyor    
        mtrec = np.array([0])
    if ttrec is None: #ttrec matrisi aynı şeyi zaman için yapıyor.
        ttrec = np.array([0])
    grad = (yf - yfprev) / dt
    #print(grad)
    
    peak = checksign(grad)
    if yk < 0.5:
        peak = 0  
            
    
    if peak == 1:
        mtrec = np.append(mtrec, yk)
        ttrec = np.append(ttrec, t)
        #print("gradclosezero = %.4f" % grad)
    else:
        regt = merge_close_values(ttrec, 5)
        regm = merge_close_values(mtrec, 0.2)
        last_peak_value = regm[-1]
        last_peak_time = regt[-1]
        founded_peak = len(regt) - 1
        
        #print("Last peak value : %.1f" % last_peak_value)
        #print("Last peak time: %.1f" % last_peak_time)
        #print("Founded peak: %.1f" % founded_peak)
        #print("gradfarfromzero = %.4f" % grad)

    return mtrec, ttrec, last_peak_value, last_peak_time, founded_peak


def lowpass_filter_single_value(signal_value, past_values=None):
    
    num_taps=71
    if past_values is None:
        past_values = np.zeros(num_taps - 1)  # İlk başta önceki değerler sıfır olacak
    cutoff_freq = 0.1  # Kesme frekansı (Hz)
    fs = 1000
    nyquist_freq = fs / 2
    norm_cutoff_freq = cutoff_freq / nyquist_freq
    taps = np.sinc(2 * norm_cutoff_freq * (np.arange(num_taps) - (num_taps - 1) / 2))
    taps *= np.hamming(num_taps)  # Hamming pencereleme
    taps /= np.sum(taps)
    
    # Önceki değerler ile filtreyi uygula
    past_values = np.append(past_values, signal_value)
    filtered_value = np.sum(past_values[-num_taps:] * taps)
    
    # Güncellenmiş past_values'i döndür
    return filtered_value, past_values[1:]


def customFilter(sample, x_array=None):
    b = np.full(30, 1/30)  # FIR filtresinin katsayıları
    
    
    
    if x_array is None:
        x_array = np.zeros_like(b)
    
    # x_array dizisini kaydırın
    x_array = np.roll(x_array, 1)
    x_array[-1] = sample
    
    
    filteredSignal = np.sum(x_array * b)
    
    return filteredSignal, x_array


def collectedmatris(deger, collected_matris):
    if not collected_matris:  
        collected_matris.append(deger)
    elif collected_matris[-1] != deger:  
        collected_matris.append(deger)
    return collected_matris

def update_PIDparameters(peakvaluematris, peaktimematris, umin, umax, ytunematris, utunematris,
                         method = 'zieglernicholas'):
    peakvaluematris = np.array(peakvaluematris)
    peaktimematris = np.array(peaktimematris)
    realpeaks = peakvaluematris[2:]
    realtimes = peaktimematris[2:]
    maxpeaks = []
    minpeaks = []
    maxpeaktimes = []
    minpeaktimes = []
    for i in range(len(realpeaks)):
        if i % 2 == 0:
            maxpeaks.append(realpeaks[i])
            maxpeaktimes.append(realtimes[i])
        else:
            minpeaks.append(realpeaks[i])
            minpeaktimes.append(realtimes[i])
    #print(minpeaks)      
    even = min(maxpeaks)
    odd = min(minpeaks)
    _, _ = np.sort([odd, even])
    ymax = max(ytunematris)
    ll = len(ytunematris)
    ymin = min(ytunematris[int(ll/4):(int(ll/4)*3)])
    print(ymin, ymax)
    
    PU = (np.mean(np.diff(maxpeaktimes)) + np.mean(np.diff(minpeaktimes))) / 2
    KU = (2 * (umax + umin) ) / (np.pi * (ymax - ymin) / 2)
    kp = (trapzoidat(ytunematris, 0.5)) / (trapzoidat(utunematris, 0.5)) # 0.5 = dt
    amp = (ymax - ymin) 
    relayheight = umax - umin
    tau = (PU / 2) / np.log((kp*relayheight + amp) / (kp*relayheight - amp))
    delay = ((PU / 2) * np.log((kp*relayheight + amp) / (kp*relayheight))) / (np.log((kp*relayheight + amp) / (kp*relayheight - amp)))   
    
    #TF ile hesaplar gelecek.. yukarıya dt değerini farklı verebilirim... 0.5=dt
    if method == 'zieglernicholas':
        KC = KU / 1.7
        TI = PU / 2.0
        TD = PU / 8.0
        
    elif method == 'tyreusluyben':
        KC = KU / 2.2
        TI = PU * 2.0
        TD = PU / 6.3
        
    elif method == 'cohencoon':
        KC = (tau / (kp * delay)) * ((16 + (3 * delay / tau)) / 12)
        TI = (delay * (32 + (6 * delay / tau))) / (13 + (8 * delay / tau))
        TD = (4 * delay) / (11 + (2 * delay / tau))
        
    elif method == 'itaedist':
        KC = (1.357 * (delay / tau) ** (-0.947)) / kp
        TI = (tau) / (0.842 * (delay / tau) ** (-0.738))
        TD = tau * 0.381 * (delay / tau) ** 0.995
        
      
    tuneMode = False     
    print("Kp founded %.2f\nTau founded %.2f\ndelay founded %.2f\nKU founded %.2f\nPU founded %.2f\n" % (kp, tau, delay, KU, PU))
    return KC, TI, TD, tuneMode
    
def relayoutput(rk, yk, umin, umax, bias, outputprev):
    s = 0
    uback = 0
    PID_back = 0
    error = rk - yk
    if yk <= (rk - bias):
        output = umax
    elif yk >= (rk + bias):
        output = umin
    else:
        output = outputprev
    
    return output, error, s, uback, PID_back 

def trapzoidat(y, dt):
    n = len(y)
    s = 0
    for k in range(n-1):
        s = s + (y[k] + y[k+1]) * (dt / 2);
        
    return s


def generateNoise(minval, maxval):
     return np.random.uniform(minval, maxval)
 

    
    
        
def checksign(value):
    # İşaret değişimini kontrol etmek için önceki değeri sakla
    if not hasattr(checksign, 'prev_value'):
        checksign.prev_value = value
        return 0
    
    # Önceki işaret ile mevcut işareti karşılaştır
    prev_sign = checksign.prev_value >= 0
    current_sign = value >= 0
    checksign.prev_value = value
    
    # İşaret değişimi varsa 1 döndür, aksi halde 0 döndür
    if current_sign != prev_sign:
        return 1
    return 0


        
class Plotter:
    def __init__(self, xlabel='X Ekseni', ylabel='Y Ekseni', title=None):
        self.x_data = []
        self.y_data = []
        self.labels = []
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.title = title
        plt.figure(figsize=(10,4))

    def addplot(self, x, y, label=None):
        """
        Verilen x ve y verilerini, belirtilen etiketle birlikte grafik üzerine çizer.
        
        Parametreler:
        x (list): x eksenindeki değerlerin listesi.
        y (list): y eksenindeki değerlerin listesi.
        label (str, optional): Veri setinin etiketi. Varsayılan olarak None.
        
        Dönüş:
        None
        """
        self.x_data.append(x)
        self.y_data.append(y)
        self.labels.append(label)

    def show(self):
        """
        Önceden eklenen tüm veri setlerini tek bir grafik üzerine çizer.
        
        Parametreler:
        None
        
        Dönüş:
        None
        """
        for i in range(len(self.x_data)):
            plt.plot(self.x_data[i], self.y_data[i], label=self.labels[i])
        plt.xlabel(self.xlabel, fontsize=14)
        plt.ylabel(self.ylabel, fontsize=14)
        plt.title(self.title, fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.grid(visible=True)
        plt.legend(fontsize=12)
        plt.show()
        
        
        
class Tablo:
    def __init__(self, *headers):
        self.headers = headers
        self.data = []

    def adddata(self, *args):
        if len(args) != len(self.headers):
            raise ValueError("Number of values provided doesn't match number of headers")
        self.data.append(args)

    def show(self):
        print(tabulate(self.data, headers=self.headers, tablefmt="presto"))

    def createpdf(self, filename='tablo.pdf'):
        pdf = PDF()
        pdf.add_page()
        pdf.set_font("Times", size=10)

        data = [self.headers] + self.data
        pdf.create_table(table_data=data, title='Error Values', cell_width='uneven', x_start=25)
        pdf.ln()
        pdf.output(filename)        




class Subplotter:
    def __init__(self, sumplot):
        self.sumplot = sumplot
        self.subplot_count = 0
        self.fig = plt.figure()


        

    def addsubplot(self, x_data, y_data, title='', xlabel='', ylabel=''):
        """
        Verilen x ve y verilerini yeni bir subplot olarak ekler.

        Parametreler:
        x_data (list): x eksenindeki değerlerin listesi.
        y_data (list): y eksenindeki değerlerin listesi.
        title (str, optional): Subplot başlığı. Varsayılan olarak boştur.
        xlabel (str, optional): x ekseninin etiketi. Varsayılan olarak boştur.
        ylabel (str, optional): y ekseninin etiketi. Varsayılan olarak boştur.

        Dönüş:
        None
        """
        self.subplot_count += 1
        ax = self.fig.add_subplot(self.sumplot, 1, self.subplot_count)
        for i in range(len(x_data)):
            ax.plot(x_data[i], y_data[i])
        ax.set_title(title, fontsize=14)
        ax.set_xlabel(xlabel, fontsize=14)
        ax.set_ylabel(ylabel, fontsize=14)        
        #ax.tick_params(axis='both', which='major', labelsize=12)
        
        ax.grid(visible=True)
       

    def show(self):
        """
        Eklenen tüm subplotları tek bir sütun halinde gösterir.

        Parametreler:
        None

        Dönüş:
        None
        """
        plt.tight_layout()
        plt.show()