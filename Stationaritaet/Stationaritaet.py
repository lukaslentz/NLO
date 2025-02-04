import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class Harmonic_Oscillator():
    def __init__(self, D=0.05, omega0=1):
        self.D = D
        self.omega0 = omega0
        self.omegad = omega0*np.sqrt(1-D**2)
        self.delta = 2 * D * omega0

    def get_t_trans(self, eps=1e-14):
        return -np.log(eps)/(self.D*self.omega0)

    def get_particular_solution(self, Omega, f):
        c_p = f / np.sqrt((self.omega0**2 - Omega**2)**2 + (self.delta * Omega)**2) 
        phi_p = np.arctan2(-self.delta * Omega , self.omega0**2 - Omega**2)
        return (
            lambda t: c_p * np.cos(Omega * t + phi_p), 
            lambda t: -c_p * Omega * np.sin(Omega * t + phi_p)
        )

    def combine_particular_solutions(self, Omegas, fs):
        Omegas = self.ensure_list(Omegas)  
        fs = self.ensure_list(fs) 
        return (
            lambda t: sum(self.get_particular_solution(Omega, f)[0](t) for Omega, f in zip(Omegas, fs)),  
            lambda t: sum(self.get_particular_solution(Omega, f)[1](t) for Omega, f in zip(Omegas, fs))  
        )

    def get_homogenous_solution(self, y0):
        c_h = np.sqrt(y0[0]**2 + ((y0[1] + self.D * self.omega0 * y0[0]) / self.omegad)**2)
        phi_h = np.arctan2(-(y0[1] + self.D * self.omega0 * y0[0]) / self.omegad, y0[0])
        return lambda t: c_h * np.exp(-self.D * self.omega0 * t) * np.cos(self.omegad * t + phi_h)

    def get_full_solution(self, Omegas, fs, y0):
        Omegas = self.ensure_list(Omegas)  
        fs = self.ensure_list(fs) 
        x_p, v_p = self.combine_particular_solutions(Omegas, fs)
        y0_h = (y0[0] - x_p(0), y0[1] - v_p(0))
        x_h = self.get_homogenous_solution(y0_h)
        return lambda t: x_p(t) + x_h(t)
    
    def rhs(self, t, state, Omegas, fs):
        Omegas = self.ensure_list(Omegas)  
        fs = self.ensure_list(fs)
        return np.array([
            state[1],
            - self.delta * state[1] 
            - self.omega0**2 * state[0] 
            + sum(f * np.cos(Omega * t) for Omega, f in zip(Omegas, fs))
        ])
    
    def create_case(self, Omegas, fs, y0):
        return(
            self.get_full_solution(Omegas, fs, y0),
            lambda t,state: self.rhs(t, state, Omegas, fs)
        )
    
    @staticmethod
    def ensure_list(value):
        if isinstance(value, (int, float)):
            return [value]
        return value
    
class Solver():
    def __init__(self, rhs, ct, digits=10):
        self.rhs = rhs
        self.ct = ct
        self.digits = digits

        self.iter = 0
        self.stationary = False
        self.points = []
        self.times = []
        self.pattern = []
        self.members = dict()
        self.member_cnt = 0

        self.stationary_pattern = []
        self.stationary_pattern_points = []
        self.stationary_pattern_ind = None
        self.pattern_inds = None
        self.stationary_pattern_times = []
        self.stationary_time = None

        self.max_sorted_pattern = []
        self.max_sorted_pattern_points = []
        self.max_sorted_pattern_times = []

        self.inits = []

        self.pattern_reps = 3 

        self.event = lambda t, state: state[1]
        #self.event.direction = -1
    
    def advance_solution(self, y0):
        t0 = self.iter * self.ct
        self.iter += 1
        te = self.iter * self.ct
        return self.integrate(t_span=[t0, te], y0=y0, t_eval=np.linspace(t0, te, 99), events=self.event)

    def integrate(self, t_span = [0, 1], y0=(1, 0), t_eval=np.linspace(0, 1, 99), events=None):
        return solve_ivp(
            self.rhs,
            t_span = t_span,
            y0 = y0,
            t_eval = t_eval,
            events = events,
            rtol = 1e-13,
            atol = 1e-13
            )
    
    def get_transient_time(self, y0=(1,0), max_iter=10):
        print('Start to find stationary pattern')
        print('-'*80)
        while not self.stationary and self.iter < max_iter:
            print(f"\t this is iteration {self.iter + 1:>6}")
            self.sol = self.advance_solution(y0)
            if self.sol.success:
                if self.sol.t_events[0].size > 0:
                    self.process_event_data(self.sol)
                y0 = self.sol.y[:,-1]
        if self.stationary:
            self.inits = [y0, self.sol.t[-1]]
            print('-'*80)
            print("successfully finished")

    def process_event_data(self, sol):
        new_points = sol.y_events[0][:,0].tolist()
        new_times = sol.t_events[0].tolist()
        pattern = self.process_points(new_points)
        self.is_subpattern(pattern)
        self.pattern.extend(pattern)

        self.points.extend(new_points)
        self.times.extend(new_times)

        self.new_points = new_points
        self.new_times = new_times

    def process_points(self, points):
        for member in {round(point, self.digits) for point in points}:
            if not member in self.members:
                self.member_cnt += 1
                self.members.setdefault(member, self.member_cnt)
        return [self.members[round(point, self.digits)] for point in points]
    
    def is_subpattern(self, test_pattern):
        pattern_inds = []
        n = len(self.pattern)
        m = len(test_pattern)

        if m > n: return []

        for i in range(n - m + 1): 
            if self.pattern[i : i + m] == test_pattern:  
                pattern_inds.append(i) 
                if len(pattern_inds) >= 2 and abs(self.times[pattern_inds[-1]] - self.times[pattern_inds[-2]] - self.ct) > 1e-9:
                    pattern_inds = []

        if len(pattern_inds) < self.pattern_reps + 1: return[] 

        eval_inds = pattern_inds[ -(self.pattern_reps + 1) : ] 
        if not self.has_equal_differences(eval_inds): return[]  

        eval_patterns = [self.pattern[ind:eval_inds[i+1]] for i,ind in enumerate(eval_inds[:-2])]
        if all(eval_patterns[i] == eval_patterns[i + 1] for i in range(len(eval_patterns) - 1) ):
            self.stationary = True
            self.stationary_pattern = eval_patterns[0]
            self.stationary_pattern_times = self.times[eval_inds[0]:eval_inds[1]]
            self.stationary_pattern_ind = eval_inds[0]
            self.pattern_inds = eval_inds
            self.process_solution()
    
    def process_solution(self):
        print("\nstart to evaluate results")
        print("\tshift pattern to the smallest possible time")
        shifted_pattern = self.stationary_pattern
        shifted_times = self.stationary_pattern_times
        shifted_ind = self.stationary_pattern_ind
        for _ in range(len(self.stationary_pattern)):
            if shifted_ind - 1 >= 0 and self.pattern[shifted_ind - 1] == shifted_pattern[-1] and abs(self.times[shifted_ind - 1] - shifted_times[-1] - self.ct) <= 1e-9:
                del shifted_pattern[-1]
                del shifted_times[-1]
                shifted_pattern.insert(0,self.pattern[shifted_ind-1])
                shifted_times.insert(0,self.times[shifted_ind-1])
                shifted_ind -= 1
            else: break
        print(f"\tshifted pattern {self.stationary_pattern_ind-shifted_ind} positions")
        self.stationary_time = shifted_times[0]
        print(f"\tsystem is stationary after {self.stationary_time} seconds")
        self.stationary_pattern_times = shifted_times
        self.stationary_pattern = shifted_pattern
        self.stationary_pattern_points = self.points[shifted_ind : shifted_ind + len(shifted_pattern)]
        print(f"\tthe first stationary pattern is")
        for point in self.stationary_pattern_points:
            print(f"\t\t{point:>9.6f}")
        max_ind = np.argmax(self.stationary_pattern_points)
        self.max_sorted_pattern = self.pattern[shifted_ind + max_ind : shifted_ind + max_ind + len(shifted_pattern)]
        self.max_sorted_pattern_points = self.points[shifted_ind + max_ind : shifted_ind + max_ind + len(shifted_pattern)]
        self.max_sorted_pattern_times = self.times[shifted_ind + max_ind :shifted_ind + max_ind + len(shifted_pattern)]

    def get_pattern_points(self, ind):
        return(
            self.points[self.pattern_inds[ind]:self.pattern_inds[ind + 1]],
            self.times[self.pattern_inds[ind]:self.pattern_inds[ind + 1]]
        )
    
    @staticmethod
    def has_equal_differences(arr, eps=1e-3):
        if not isinstance(arr, np.ndarray): arr = np.array(arr) 
        differences = np.diff(arr)
        return np.all(np.abs(differences - differences[0]) < eps)

def plot_line(ax, x, y, color='blue',linestyle='-'):
    ax.plot(x, y, color=color, linestyle=linestyle, zorder=5)
    # ax.grid(True)
    return ax

def plot_points(ax, x, y, color='blue'):
    ax.scatter(x, y, color=color, zorder=10)
    # ax.grid(True)
    return ax

def main():

    # Testcase definieren
    harmonic_oscillator = Harmonic_Oscillator(D=0.01, omega0 = 1.65)
    Omegas = [1, 1/2] # Anzahl der Omegas legt fest, wieviele Frequenzen in der Antwort enthalten sind
    fs = (1, 1) # Amplituden der Anregungen mit Erregerkreisfrequenzen Omegas
    y0 = (10, 2)

    (x, rhs) = harmonic_oscillator.create_case(Omegas, fs, y0)
    
    print(f"\ntheoretical limit for transient time {harmonic_oscillator.get_t_trans(eps=1e-12)} s\n")
  
    # Einschwingen berechnen
    ct =  1 * 2 * np.pi / min(Omegas) # Zeit mit der abgetatstet wird
    solver = Solver(rhs, ct, digits = 7)
    solver.get_transient_time(y0 = y0, max_iter = 1e10)

    te = solver.times[-1] # letzte Zeit bis zu der gerechnet wurde 
    ts = solver.stationary_time # Zeit ab der das System stationär wird

    # Graphische Auswertung
    times = np.linspace(0, te, 9999)
    
    _, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.grid(True, zorder=0)
    ax.set_xlim(-0.1*ct, te+.5*ct)
    ax.set_xticks(np.arange(0, te+ct, ct))

    plot_line(ax, times, x(times)) # Verlauf des Signals plotten

    xp, t = solver.points, solver.times
    plot_points(ax, t, xp, color='black')

    colors = ['brown', 'orange', 'purple', 'orange','green','black', 'brown', 'orange', 'purple', 'orange','green']

    for i in range(0, len(solver.pattern_inds) -1):
        pattern_x, pattern_t = solver.get_pattern_points(i)
        plot_points(ax, pattern_t, pattern_x, color=colors[i]) # Gefundene, sich wiederhohlende Muster plotten

    #plot_points(ax, solver.stationary_pattern_times, solver.stationary_pattern_points, color='pink') # Die frühste variante der Muster plotten

    ax.text(
        0.3, 1.08,
        f"Zeit bis zur Stationarität beträgt: {ts:>8.2f} s \n und die maximale stationäre Amplitude ist {max(solver.stationary_pattern_points)}",
        fontsize=12, 
        ha='left', 
        va='center',  
        transform=ax.transAxes,
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')
        )
    
    ax.axvline(x = ts, color='red', linestyle='--', linewidth=2) # stationäre Zeit einzeichen
    
    # Ränder des letzten Bereichs der berechnet wurde
    ax.axvline(x = solver.new_times[0], color='green', linestyle='--', linewidth=2)
    ax.axvline(x = solver.new_times[-1], color='green', linestyle='--', linewidth=2)

    plot_points(ax, solver.new_times, solver.new_points, color='red') # Wendepunkte im letzten Bereich, die zum finden des Musters verwendet wurden

    plt.show()



if __name__ == "__main__":
    main()
