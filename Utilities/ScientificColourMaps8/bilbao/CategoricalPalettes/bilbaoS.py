# 
#         bilbaoS
#                   www.fabiocrameri.ch/colourmaps
from matplotlib.colors import LinearSegmentedColormap      
      
cm_data = [[0.29914, 0.00032437, 0.0030268],      
           [0.97154, 0.97151, 0.97136],      
           [0.66084, 0.50434, 0.36855],      
           [0.59667, 0.30066, 0.31369],      
           [0.74569, 0.71659, 0.57883],      
           [0.46742, 0.15698, 0.17771],      
           [0.63663, 0.41478, 0.34915],      
           [0.68821, 0.60091, 0.39598],      
           [0.80036, 0.79399, 0.7641],      
           [0.6223, 0.36343, 0.33749],      
           [0.54063, 0.22756, 0.25552],      
           [0.71532, 0.66509, 0.46623],      
           [0.3915, 0.087894, 0.10233],      
           [0.86624, 0.86518, 0.86021],      
           [0.67276, 0.54821, 0.37816],      
           [0.77093, 0.75498, 0.67909],      
           [0.64929, 0.46161, 0.35929],      
           [0.62987, 0.38996, 0.3437],      
           [0.67972, 0.57323, 0.38456],      
           [0.91647, 0.91623, 0.91512],      
           [0.50497, 0.19208, 0.21755],      
           [0.34614, 0.044376, 0.062581],      
           [0.64304, 0.43852, 0.35429],      
           [0.70001, 0.63245, 0.42126],      
           [0.75871, 0.73659, 0.63108],      
           [0.61235, 0.33406, 0.32876],      
           [0.57236, 0.26419, 0.28889],      
           [0.73123, 0.69345, 0.52255],      
           [0.66631, 0.5245, 0.37292],      
           [0.42814, 0.12109, 0.13725],      
           [0.82527, 0.82207, 0.8071],      
           [0.65468, 0.48155, 0.36361],      
           [0.7842, 0.77372, 0.72423],      
           [0.6333, 0.40254, 0.34648],      
           [0.68366, 0.58662, 0.38919],      
           [0.77722, 0.76409, 0.70185],      
           [0.70737, 0.64901, 0.4415],      
           [0.62624, 0.37698, 0.34076],      
           [0.6695, 0.53624, 0.37548],      
           [0.55721, 0.24572, 0.27303],      
           [0.73867, 0.70555, 0.5511],      
           [0.60547, 0.31788, 0.3223],      
           [0.36925, 0.067261, 0.083028],      
           [0.6178, 0.34919, 0.33364],      
           [0.58563, 0.28265, 0.30255],      
           [0.32267, 0.020224, 0.035143],      
           [0.52312, 0.2097, 0.23691],      
           [0.89054, 0.89002, 0.88753],      
           [0.72339, 0.68, 0.49384],      
           [0.84424, 0.84231, 0.83327],      
           [0.63987, 0.42676, 0.35174],      
           [0.94356, 0.94347, 0.94302],      
           [0.48638, 0.17457, 0.19775],      
           [0.76486, 0.7459, 0.65557],      
           [0.65776, 0.49292, 0.36608],      
           [0.64618, 0.45012, 0.35682],      
           [0.75235, 0.72687, 0.60551],      
           [0.67613, 0.56049, 0.38109],      
           [0.69359, 0.61623, 0.40626],      
           [0.44804, 0.13919, 0.15745],      
           [0.40753, 0.10245, 0.11715],      
           [0.66318, 0.51293, 0.3704],      
           [0.80976, 0.80488, 0.78206],      
           [0.7926, 0.78457, 0.74676],      
           [0.65237, 0.47301, 0.36177],      
           [0.95745, 0.9574, 0.95713],      
           [0.74223, 0.71118, 0.56508],      
           [0.53197, 0.2186, 0.24631],      
           [0.6593, 0.49862, 0.3673],      
           [0.69077, 0.60845, 0.4006],      
           [0.78058, 0.7688, 0.71306],      
           [0.67442, 0.5543, 0.37958],      
           [0.92989, 0.92974, 0.92902],      
           [0.45779, 0.14812, 0.16759],      
           [0.60913, 0.32609, 0.32575],      
           [0.64773, 0.45588, 0.35806],      
           [0.47695, 0.16576, 0.18778],      
           [0.57926, 0.27344, 0.29602],      
           [0.54906, 0.23661, 0.26444],      
           [0.7113, 0.65717, 0.45339],      
           [0.62809, 0.38354, 0.34227],      
           [0.75556, 0.73179, 0.61844],      
           [0.6679, 0.53035, 0.3742],      
           [0.59143, 0.29174, 0.30845],      
           [0.71937, 0.67273, 0.4798],      
           [0.33444, 0.031738, 0.050036],      
           [0.6779, 0.56679, 0.38273],      
           [0.90334, 0.90298, 0.90129],      
           [0.63497, 0.4087, 0.34783],      
           [0.67111, 0.54219, 0.37681],      
           [0.83434, 0.83182, 0.82006],      
           [0.60133, 0.30941, 0.31829],      
           [0.76788, 0.75046, 0.66744],      
           [0.49572, 0.18331, 0.2077],      
           [0.78816, 0.77894, 0.73544],      
           [0.62015, 0.35642, 0.33566],      
           [0.7618, 0.74129, 0.64345],      
           [0.73501, 0.69965, 0.5369],      
           [0.43817, 0.13026, 0.14735],      
           [0.87815, 0.87739, 0.87383]]      
      
bilbaoS_map = LinearSegmentedColormap.from_list('bilbaoS', cm_data)      
# For use of "viscm view"      
test_cm = bilbaoS_map      
      
if __name__ == "__main__":      
    import matplotlib.pyplot as plt      
    import numpy as np      
      
    try:      
        from viscm import viscm      
        viscm(bilbaoS_map)      
    except ImportError:      
        print("viscm not found, falling back on simple display")      
        plt.imshow(np.linspace(0, 100, 256)[None, :], aspect='auto',      
                   cmap=bilbaoS_map)      
    plt.show()      