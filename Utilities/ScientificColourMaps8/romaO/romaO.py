# 
#         romaO
#                   www.fabiocrameri.ch/colourmaps
from matplotlib.colors import LinearSegmentedColormap      
      
cm_data = [[0.45137, 0.22346, 0.34187],      
           [0.45418, 0.22244, 0.3361],      
           [0.45696, 0.22158, 0.33043],      
           [0.45975, 0.2209, 0.32483],      
           [0.46251, 0.22035, 0.31935],      
           [0.46527, 0.21994, 0.31394],      
           [0.46803, 0.21968, 0.30862],      
           [0.47078, 0.21958, 0.30337],      
           [0.47352, 0.21962, 0.29822],      
           [0.47628, 0.21982, 0.29316],      
           [0.47902, 0.22017, 0.28818],      
           [0.48178, 0.22067, 0.2833],      
           [0.48453, 0.2213, 0.2785],      
           [0.48731, 0.22208, 0.27379],      
           [0.49008, 0.22304, 0.26917],      
           [0.49286, 0.22411, 0.26461],      
           [0.49567, 0.22536, 0.26016],      
           [0.4985, 0.22677, 0.25579],      
           [0.50134, 0.22833, 0.25153],      
           [0.50419, 0.22999, 0.24733],      
           [0.50707, 0.23188, 0.24322],      
           [0.50997, 0.23387, 0.23923],      
           [0.5129, 0.23605, 0.23533],      
           [0.51584, 0.23835, 0.23151],      
           [0.51884, 0.24082, 0.22779],      
           [0.52184, 0.24345, 0.22414],      
           [0.52489, 0.24625, 0.22065],      
           [0.52797, 0.2492, 0.2172],      
           [0.53108, 0.25231, 0.21387],      
           [0.53423, 0.25556, 0.21064],      
           [0.53742, 0.25899, 0.20753],      
           [0.54063, 0.26255, 0.20452],      
           [0.54389, 0.26628, 0.20158],      
           [0.54718, 0.27017, 0.19879],      
           [0.55051, 0.27419, 0.19613],      
           [0.55389, 0.27839, 0.19356],      
           [0.55731, 0.28273, 0.19109],      
           [0.56075, 0.2872, 0.18877],      
           [0.56424, 0.29186, 0.18655],      
           [0.56777, 0.29665, 0.18446],      
           [0.57134, 0.30157, 0.18248],      
           [0.57495, 0.30666, 0.18065],      
           [0.5786, 0.31186, 0.17898],      
           [0.58228, 0.31724, 0.17743],      
           [0.58602, 0.32275, 0.17597],      
           [0.58977, 0.32838, 0.17473],      
           [0.59358, 0.33415, 0.17358],      
           [0.59742, 0.34005, 0.17261],      
           [0.60129, 0.34606, 0.17179],      
           [0.60519, 0.35223, 0.17114],      
           [0.60915, 0.35851, 0.17065],      
           [0.61311, 0.36491, 0.17034],      
           [0.61713, 0.37143, 0.1702],      
           [0.62118, 0.37808, 0.17023],      
           [0.62526, 0.38483, 0.17046],      
           [0.62937, 0.39171, 0.17087],      
           [0.63352, 0.39869, 0.17148],      
           [0.63769, 0.40579, 0.17229],      
           [0.6419, 0.41299, 0.17332],      
           [0.64613, 0.42029, 0.17458],      
           [0.65041, 0.42771, 0.176],      
           [0.6547, 0.43522, 0.17774],      
           [0.65904, 0.44283, 0.17962],      
           [0.66341, 0.45054, 0.18175],      
           [0.6678, 0.45834, 0.18416],      
           [0.67222, 0.46625, 0.1868],      
           [0.67667, 0.47425, 0.18968],      
           [0.68114, 0.48233, 0.19283],      
           [0.68566, 0.49051, 0.19624],      
           [0.69019, 0.49878, 0.19987],      
           [0.69474, 0.50712, 0.20384],      
           [0.69933, 0.51554, 0.20803],      
           [0.70394, 0.52406, 0.21251],      
           [0.70858, 0.53265, 0.21726],      
           [0.71322, 0.5413, 0.22229],      
           [0.7179, 0.55003, 0.22761],      
           [0.72257, 0.55881, 0.23318],      
           [0.72727, 0.56767, 0.23907],      
           [0.73197, 0.57658, 0.24521],      
           [0.73666, 0.58553, 0.25168],      
           [0.74136, 0.59451, 0.25837],      
           [0.74605, 0.60354, 0.26537],      
           [0.75073, 0.61259, 0.27263],      
           [0.75538, 0.62166, 0.28017],      
           [0.76001, 0.63075, 0.28796],      
           [0.7646, 0.63982, 0.29602],      
           [0.76914, 0.64889, 0.30433],      
           [0.77363, 0.65793, 0.31287],      
           [0.77806, 0.66694, 0.32165],      
           [0.78242, 0.6759, 0.33066],      
           [0.78669, 0.68481, 0.33988],      
           [0.79087, 0.69365, 0.34929],      
           [0.79494, 0.7024, 0.35888],      
           [0.7989, 0.71106, 0.36867],      
           [0.80273, 0.71961, 0.37859],      
           [0.80642, 0.72803, 0.38866],      
           [0.80996, 0.73631, 0.39885],      
           [0.81334, 0.74446, 0.40916],      
           [0.81655, 0.75244, 0.41957],      
           [0.81956, 0.76025, 0.43004],      
           [0.82239, 0.76787, 0.44057],      
           [0.82501, 0.7753, 0.45115],      
           [0.82742, 0.78252, 0.46174],      
           [0.8296, 0.78953, 0.47235],      
           [0.83155, 0.79631, 0.48293],      
           [0.83326, 0.80287, 0.49349],      
           [0.83472, 0.80919, 0.50402],      
           [0.83592, 0.81526, 0.51449],      
           [0.83686, 0.82109, 0.52487],      
           [0.83753, 0.82666, 0.53517],      
           [0.83793, 0.83198, 0.54537],      
           [0.83805, 0.83703, 0.55546],      
           [0.83788, 0.84182, 0.56542],      
           [0.83744, 0.84635, 0.57525],      
           [0.8367, 0.85061, 0.58493],      
           [0.83567, 0.85462, 0.59446],      
           [0.83435, 0.85835, 0.60382],      
           [0.83274, 0.86183, 0.61301],      
           [0.83084, 0.86504, 0.62202],      
           [0.82864, 0.868, 0.63085],      
           [0.82615, 0.87068, 0.63949],      
           [0.82337, 0.87312, 0.64792],      
           [0.8203, 0.87531, 0.65617],      
           [0.81695, 0.87724, 0.6642],      
           [0.81331, 0.87892, 0.67203],      
           [0.80939, 0.88036, 0.67964],      
           [0.80518, 0.88156, 0.68705],      
           [0.80071, 0.8825, 0.69424],      
           [0.79595, 0.88322, 0.70121],      
           [0.79094, 0.8837, 0.70797],      
           [0.78566, 0.88395, 0.7145],      
           [0.78012, 0.88396, 0.72082],      
           [0.77433, 0.88375, 0.72692],      
           [0.7683, 0.88331, 0.73279],      
           [0.76203, 0.88264, 0.73844],      
           [0.75553, 0.88177, 0.74387],      
           [0.74879, 0.88066, 0.74908],      
           [0.74184, 0.87934, 0.75407],      
           [0.73468, 0.87781, 0.75884],      
           [0.72731, 0.87607, 0.76339],      
           [0.71976, 0.87411, 0.76772],      
           [0.71201, 0.87195, 0.77184],      
           [0.70408, 0.86958, 0.77573],      
           [0.69599, 0.86701, 0.77941],      
           [0.68774, 0.86425, 0.78288],      
           [0.67934, 0.86127, 0.78614],      
           [0.67081, 0.85811, 0.78919],      
           [0.66215, 0.85476, 0.79202],      
           [0.65336, 0.8512, 0.79465],      
           [0.64448, 0.84747, 0.79707],      
           [0.6355, 0.84356, 0.7993],      
           [0.62645, 0.83947, 0.80131],      
           [0.61732, 0.83519, 0.80313],      
           [0.60814, 0.83075, 0.80476],      
           [0.59891, 0.82614, 0.80619],      
           [0.58965, 0.82137, 0.80743],      
           [0.58037, 0.81644, 0.80848],      
           [0.57108, 0.81135, 0.80935],      
           [0.56181, 0.80612, 0.81004],      
           [0.55255, 0.80074, 0.81055],      
           [0.54332, 0.79522, 0.81088],      
           [0.53412, 0.78958, 0.81105],      
           [0.525, 0.7838, 0.81105],      
           [0.51593, 0.77791, 0.81088],      
           [0.50695, 0.77189, 0.81055],      
           [0.49808, 0.76577, 0.81007],      
           [0.48928, 0.75954, 0.80944],      
           [0.48061, 0.75321, 0.80866],      
           [0.47207, 0.7468, 0.80773],      
           [0.46365, 0.74029, 0.80667],      
           [0.45539, 0.7337, 0.80546],      
           [0.44728, 0.72703, 0.80413],      
           [0.43934, 0.7203, 0.80266],      
           [0.43158, 0.7135, 0.80107],      
           [0.42398, 0.70664, 0.79936],      
           [0.41658, 0.69971, 0.79752],      
           [0.40938, 0.69275, 0.79557],      
           [0.40237, 0.68572, 0.79351],      
           [0.3956, 0.67865, 0.79133],      
           [0.38903, 0.67155, 0.78905],      
           [0.38267, 0.66441, 0.78666],      
           [0.37656, 0.65724, 0.78416],      
           [0.37066, 0.65003, 0.78155],      
           [0.36502, 0.64279, 0.77884],      
           [0.35961, 0.63552, 0.77604],      
           [0.35446, 0.62824, 0.77312],      
           [0.34955, 0.62094, 0.77011],      
           [0.3449, 0.6136, 0.767],      
           [0.34051, 0.60625, 0.76378],      
           [0.33637, 0.59889, 0.76047],      
           [0.33253, 0.59151, 0.75704],      
           [0.32893, 0.58412, 0.75351],      
           [0.32559, 0.57671, 0.74987],      
           [0.32256, 0.56928, 0.74613],      
           [0.31978, 0.56186, 0.74228],      
           [0.31727, 0.55441, 0.7383],      
           [0.31505, 0.54695, 0.73422],      
           [0.31311, 0.53948, 0.73002],      
           [0.31144, 0.53201, 0.72569],      
           [0.31007, 0.52453, 0.72124],      
           [0.30897, 0.51704, 0.71667],      
           [0.30811, 0.50955, 0.71197],      
           [0.30755, 0.50205, 0.70713],      
           [0.30726, 0.49456, 0.70216],      
           [0.30723, 0.48707, 0.69706],      
           [0.30746, 0.47958, 0.69182],      
           [0.30795, 0.4721, 0.68643],      
           [0.3087, 0.46463, 0.6809],      
           [0.30968, 0.45716, 0.67525],      
           [0.31088, 0.44973, 0.66944],      
           [0.31228, 0.44232, 0.6635],      
           [0.31393, 0.43493, 0.65741],      
           [0.31578, 0.42758, 0.65118],      
           [0.3178, 0.42025, 0.64482],      
           [0.32001, 0.41299, 0.63833],      
           [0.32238, 0.40577, 0.6317],      
           [0.32489, 0.39861, 0.62495],      
           [0.32755, 0.39152, 0.61809],      
           [0.33035, 0.38448, 0.61111],      
           [0.33327, 0.37755, 0.60402],      
           [0.33627, 0.37068, 0.59684],      
           [0.33939, 0.36392, 0.58955],      
           [0.34257, 0.35728, 0.58219],      
           [0.3458, 0.35073, 0.57476],      
           [0.34912, 0.34428, 0.56727],      
           [0.35247, 0.33797, 0.55971],      
           [0.35587, 0.33179, 0.55212],      
           [0.35927, 0.32574, 0.54448],      
           [0.36271, 0.31986, 0.53684],      
           [0.36617, 0.31411, 0.52917],      
           [0.36961, 0.30852, 0.52148],      
           [0.37306, 0.30306, 0.51382],      
           [0.37652, 0.2978, 0.50615],      
           [0.37994, 0.29269, 0.49854],      
           [0.38336, 0.28775, 0.49094],      
           [0.38674, 0.28301, 0.48337],      
           [0.39011, 0.27842, 0.47586],      
           [0.39346, 0.27401, 0.4684],      
           [0.39677, 0.26978, 0.461],      
           [0.40006, 0.26573, 0.45366],      
           [0.40333, 0.26185, 0.4464],      
           [0.40655, 0.25815, 0.43921],      
           [0.40974, 0.25466, 0.43212],      
           [0.4129, 0.25132, 0.42509],      
           [0.41602, 0.24817, 0.41813],      
           [0.41912, 0.24515, 0.41128],      
           [0.42218, 0.24235, 0.40451],      
           [0.42522, 0.23972, 0.39784],      
           [0.42823, 0.23728, 0.39126],      
           [0.43121, 0.23498, 0.38475],      
           [0.43415, 0.23282, 0.37836],      
           [0.43708, 0.23086, 0.37204],      
           [0.43998, 0.22907, 0.36583],      
           [0.44286, 0.22743, 0.3597],      
           [0.44571, 0.22596, 0.35366],      
           [0.44855, 0.2246, 0.34773]]      
      
romaO_map = LinearSegmentedColormap.from_list('romaO', cm_data)      
# For use of "viscm view"      
test_cm = romaO_map      
      
if __name__ == "__main__":      
    import matplotlib.pyplot as plt      
    import numpy as np      
      
    try:      
        from viscm import viscm      
        viscm(romaO_map)      
    except ImportError:      
        print("viscm not found, falling back on simple display")      
        plt.imshow(np.linspace(0, 100, 256)[None, :], aspect='auto',      
                   cmap=romaO_map)      
    plt.show()      