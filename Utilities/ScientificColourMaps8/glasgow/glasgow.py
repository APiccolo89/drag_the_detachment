# 
#         glasgow
#                   www.fabiocrameri.ch/colourmaps
from matplotlib.colors import LinearSegmentedColormap      
      
cm_data = [[0.21181, 0.073933, 0.22061],      
           [0.21584, 0.074823, 0.21741],      
           [0.2198, 0.075756, 0.21424],      
           [0.22367, 0.076703, 0.21107],      
           [0.22749, 0.07767, 0.20791],      
           [0.23125, 0.078655, 0.20477],      
           [0.23496, 0.079657, 0.20161],      
           [0.23858, 0.08067, 0.19849],      
           [0.24219, 0.081714, 0.19539],      
           [0.24573, 0.082671, 0.1923],      
           [0.24925, 0.083676, 0.18923],      
           [0.25274, 0.084681, 0.18619],      
           [0.25617, 0.085604, 0.18312],      
           [0.2596, 0.086656, 0.1801],      
           [0.26298, 0.087669, 0.17711],      
           [0.26634, 0.088633, 0.17412],      
           [0.26968, 0.089618, 0.17115],      
           [0.27298, 0.0906, 0.16819],      
           [0.27629, 0.091604, 0.16524],      
           [0.27957, 0.092519, 0.16237],      
           [0.28284, 0.093477, 0.15943],      
           [0.28607, 0.094504, 0.15658],      
           [0.28934, 0.095451, 0.15367],      
           [0.29257, 0.096347, 0.15081],      
           [0.29581, 0.097347, 0.14796],      
           [0.29906, 0.09833, 0.14514],      
           [0.30229, 0.099261, 0.14228],      
           [0.30557, 0.10022, 0.1394],      
           [0.30883, 0.10118, 0.13657],      
           [0.31209, 0.10216, 0.13376],      
           [0.31541, 0.10316, 0.13093],      
           [0.31873, 0.10412, 0.12807],      
           [0.32208, 0.10509, 0.1252],      
           [0.32545, 0.10617, 0.1223],      
           [0.32889, 0.10719, 0.11944],      
           [0.33235, 0.10818, 0.11656],      
           [0.33585, 0.10931, 0.11364],      
           [0.33941, 0.11041, 0.11071],      
           [0.34301, 0.11153, 0.10771],      
           [0.34665, 0.11271, 0.10469],      
           [0.35037, 0.1139, 0.1017],      
           [0.35412, 0.11515, 0.098667],      
           [0.35795, 0.11651, 0.095593],      
           [0.36183, 0.11789, 0.092462],      
           [0.36575, 0.11935, 0.089337],      
           [0.36973, 0.12088, 0.086106],      
           [0.37376, 0.12251, 0.082919],      
           [0.37783, 0.12428, 0.079656],      
           [0.38192, 0.12619, 0.07636],      
           [0.38604, 0.12818, 0.073124],      
           [0.39017, 0.13034, 0.069771],      
           [0.3943, 0.13258, 0.066344],      
           [0.39842, 0.13503, 0.062973],      
           [0.4025, 0.13764, 0.059485],      
           [0.40654, 0.14038, 0.055873],      
           [0.41051, 0.14333, 0.052208],      
           [0.41438, 0.14644, 0.048631],      
           [0.41814, 0.14977, 0.044938],      
           [0.42179, 0.1532, 0.041226],      
           [0.4253, 0.15687, 0.037499],      
           [0.42864, 0.16066, 0.033747],      
           [0.43181, 0.16463, 0.030454],      
           [0.43479, 0.16873, 0.027345],      
           [0.43757, 0.17293, 0.024453],      
           [0.44013, 0.17727, 0.021779],      
           [0.4425, 0.18166, 0.019323],      
           [0.44465, 0.18619, 0.017081],      
           [0.44658, 0.19073, 0.015046],      
           [0.44831, 0.19533, 0.013206],      
           [0.44984, 0.19994, 0.011537],      
           [0.45119, 0.20462, 0.0099048],      
           [0.45236, 0.20926, 0.0086143],      
           [0.45334, 0.21392, 0.0074732],      
           [0.45419, 0.21855, 0.0064682],      
           [0.45489, 0.22316, 0.0055888],      
           [0.45546, 0.22775, 0.0048194],      
           [0.45592, 0.2323, 0.0041464],      
           [0.45627, 0.23683, 0.0035565],      
           [0.45655, 0.2413, 0.003037],      
           [0.45674, 0.24574, 0.0025765],      
           [0.45688, 0.25015, 0.0021644],      
           [0.45695, 0.25456, 0.0018062],      
           [0.45698, 0.25891, 0.0015128],      
           [0.45696, 0.26322, 0.0012797],      
           [0.45692, 0.26749, 0.0011036],      
           [0.45684, 0.27176, 0.00098266],      
           [0.45674, 0.27603, 0.0009164],      
           [0.45661, 0.28022, 0.0009059],      
           [0.45647, 0.28442, 0.0009539],      
           [0.45631, 0.28862, 0.0010649],      
           [0.45613, 0.29278, 0.0012453],      
           [0.45594, 0.29695, 0.0015036],      
           [0.45573, 0.30108, 0.0018503],      
           [0.45551, 0.30524, 0.0022982],      
           [0.45527, 0.30938, 0.0028627],      
           [0.45502, 0.31349, 0.0035614],      
           [0.45474, 0.31761, 0.0044146],      
           [0.45443, 0.32172, 0.0054454],      
           [0.4541, 0.32583, 0.0066794],      
           [0.45375, 0.32996, 0.0081439],      
           [0.45337, 0.33407, 0.0098688],      
           [0.45296, 0.33816, 0.012123],      
           [0.45251, 0.34228, 0.014455],      
           [0.45203, 0.34636, 0.017184],      
           [0.4515, 0.35047, 0.020328],      
           [0.45094, 0.35454, 0.023928],      
           [0.45033, 0.3586, 0.028028],      
           [0.44969, 0.36266, 0.032682],      
           [0.449, 0.3667, 0.038112],      
           [0.44825, 0.37071, 0.043789],      
           [0.44746, 0.37472, 0.049799],      
           [0.44664, 0.3787, 0.055892],      
           [0.44576, 0.38265, 0.062145],      
           [0.44486, 0.38656, 0.068579],      
           [0.4439, 0.39046, 0.075044],      
           [0.44291, 0.39431, 0.081754],      
           [0.44187, 0.39814, 0.08848],      
           [0.44079, 0.40191, 0.095363],      
           [0.43969, 0.40568, 0.1023],      
           [0.43856, 0.40938, 0.10939],      
           [0.43741, 0.41306, 0.1165],      
           [0.43622, 0.41669, 0.12365],      
           [0.43502, 0.42027, 0.13095],      
           [0.43378, 0.42383, 0.13827],      
           [0.43255, 0.42736, 0.14564],      
           [0.4313, 0.43084, 0.153],      
           [0.43003, 0.43428, 0.16043],      
           [0.42875, 0.4377, 0.16791],      
           [0.42748, 0.44107, 0.17541],      
           [0.42618, 0.44443, 0.18291],      
           [0.4249, 0.44774, 0.19045],      
           [0.4236, 0.45105, 0.19799],      
           [0.42231, 0.45432, 0.20554],      
           [0.42102, 0.45756, 0.21309],      
           [0.41973, 0.46078, 0.22069],      
           [0.41843, 0.46399, 0.22824],      
           [0.41715, 0.46718, 0.2358],      
           [0.41586, 0.47034, 0.24334],      
           [0.41457, 0.47348, 0.25092],      
           [0.41329, 0.47663, 0.25848],      
           [0.41201, 0.47975, 0.26605],      
           [0.41073, 0.48286, 0.27361],      
           [0.40945, 0.48596, 0.28114],      
           [0.40817, 0.48904, 0.2887],      
           [0.40689, 0.49213, 0.29625],      
           [0.40563, 0.4952, 0.30378],      
           [0.40434, 0.49827, 0.31133],      
           [0.40308, 0.50133, 0.31887],      
           [0.4018, 0.50438, 0.3264],      
           [0.40055, 0.50742, 0.33395],      
           [0.39926, 0.51046, 0.34148],      
           [0.398, 0.5135, 0.34901],      
           [0.39673, 0.51653, 0.35654],      
           [0.39547, 0.51956, 0.36405],      
           [0.3942, 0.52258, 0.37156],      
           [0.39294, 0.52561, 0.37908],      
           [0.39169, 0.52863, 0.38659],      
           [0.39044, 0.53164, 0.39409],      
           [0.3892, 0.53464, 0.40159],      
           [0.38796, 0.53766, 0.40908],      
           [0.38673, 0.54066, 0.41656],      
           [0.38552, 0.54366, 0.42403],      
           [0.38434, 0.54666, 0.43151],      
           [0.38318, 0.54965, 0.43895],      
           [0.38204, 0.55265, 0.44641],      
           [0.38095, 0.55564, 0.45385],      
           [0.37989, 0.55864, 0.46129],      
           [0.37889, 0.56166, 0.46873],      
           [0.37796, 0.56466, 0.47616],      
           [0.3771, 0.56768, 0.48358],      
           [0.37634, 0.5707, 0.49102],      
           [0.37567, 0.57374, 0.49844],      
           [0.37513, 0.57679, 0.50586],      
           [0.37473, 0.57986, 0.5133],      
           [0.3745, 0.58293, 0.52074],      
           [0.37444, 0.58605, 0.52818],      
           [0.37458, 0.58916, 0.53562],      
           [0.37495, 0.59233, 0.54308],      
           [0.37556, 0.59549, 0.55053],      
           [0.37646, 0.59869, 0.55798],      
           [0.37763, 0.60192, 0.56542],      
           [0.37911, 0.60517, 0.57287],      
           [0.38093, 0.60845, 0.58029],      
           [0.38308, 0.61174, 0.58769],      
           [0.38557, 0.61506, 0.59507],      
           [0.38844, 0.6184, 0.6024],      
           [0.39167, 0.62174, 0.6097],      
           [0.39525, 0.62511, 0.61695],      
           [0.39918, 0.62848, 0.62413],      
           [0.40347, 0.63185, 0.63125],      
           [0.40808, 0.63522, 0.63828],      
           [0.41302, 0.63858, 0.64524],      
           [0.41824, 0.64194, 0.65212],      
           [0.42374, 0.64527, 0.6589],      
           [0.4295, 0.6486, 0.6656],      
           [0.43548, 0.6519, 0.6722],      
           [0.44166, 0.65518, 0.6787],      
           [0.448, 0.65844, 0.6851],      
           [0.45452, 0.66167, 0.69141],      
           [0.46115, 0.66486, 0.69762],      
           [0.46789, 0.66804, 0.70375],      
           [0.47473, 0.67118, 0.7098],      
           [0.48162, 0.6743, 0.71576],      
           [0.48858, 0.67739, 0.72163],      
           [0.49556, 0.68045, 0.72744],      
           [0.50258, 0.68348, 0.73318],      
           [0.50962, 0.68649, 0.73886],      
           [0.51666, 0.68948, 0.74447],      
           [0.5237, 0.69245, 0.75003],      
           [0.53073, 0.69538, 0.75554],      
           [0.53775, 0.6983, 0.761],      
           [0.54475, 0.7012, 0.76641],      
           [0.55172, 0.70407, 0.77179],      
           [0.55867, 0.70693, 0.77712],      
           [0.56559, 0.70977, 0.78241],      
           [0.57249, 0.71259, 0.78766],      
           [0.57935, 0.7154, 0.79289],      
           [0.58617, 0.71818, 0.79807],      
           [0.59296, 0.72094, 0.80322],      
           [0.59971, 0.72369, 0.80834],      
           [0.60643, 0.72642, 0.81341],      
           [0.6131, 0.72914, 0.81847],      
           [0.61976, 0.73184, 0.82349],      
           [0.62636, 0.73453, 0.82848],      
           [0.63295, 0.7372, 0.83345],      
           [0.6395, 0.73986, 0.83839],      
           [0.64603, 0.74251, 0.84332],      
           [0.65254, 0.74515, 0.84822],      
           [0.65903, 0.74779, 0.85312],      
           [0.66553, 0.75042, 0.858],      
           [0.67202, 0.75304, 0.86289],      
           [0.67852, 0.75568, 0.86778],      
           [0.68503, 0.75832, 0.87268],      
           [0.69157, 0.76097, 0.87759],      
           [0.69814, 0.76364, 0.88252],      
           [0.70476, 0.76632, 0.88748],      
           [0.71144, 0.76902, 0.89249],      
           [0.71818, 0.77176, 0.89753],      
           [0.72499, 0.77451, 0.90262],      
           [0.73189, 0.77731, 0.90777],      
           [0.73888, 0.78014, 0.91297],      
           [0.74597, 0.78301, 0.91825],      
           [0.75317, 0.78593, 0.9236],      
           [0.7605, 0.7889, 0.92902],      
           [0.76795, 0.79191, 0.93452],      
           [0.77554, 0.79498, 0.94011],      
           [0.78326, 0.7981, 0.94578],      
           [0.79113, 0.80128, 0.95154],      
           [0.79915, 0.80452, 0.95738],      
           [0.8073, 0.80782, 0.9633],      
           [0.8156, 0.81117, 0.96929],      
           [0.82403, 0.81457, 0.97536],      
           [0.8326, 0.81802, 0.9815],      
           [0.84128, 0.82153, 0.9877],      
           [0.85008, 0.82507, 0.99395],      
           [0.85897, 0.82865, 1]]      
      
glasgow_map = LinearSegmentedColormap.from_list('glasgow', cm_data)      
# For use of "viscm view"      
test_cm = glasgow_map      
      
if __name__ == "__main__":      
    import matplotlib.pyplot as plt      
    import numpy as np      
      
    try:      
        from viscm import viscm      
        viscm(glasgow_map)      
    except ImportError:      
        print("viscm not found, falling back on simple display")      
        plt.imshow(np.linspace(0, 100, 256)[None, :], aspect='auto',      
                   cmap=glasgow_map)      
    plt.show()      
