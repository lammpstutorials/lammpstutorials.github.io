import numpy as np

colors = {
  "myblue": [0.0, 0.588, 0.694],
  "mydarkblue": [0.008, 0.294, 0.478],
  "mycyan": [0.267, 0.718, 0.761],
  "myorange": [1.0, 0.682, 0.286],
  "mylightgray": [0.1, 0.1, 0.1],
  "mygray": [0.5, 0.5, 0.5],
  "mydarkgray": [0.9, 0.9, 0.9],
  "myblueANR": [0.224, 0.467, 0.878],
  "myredANR": [0.867, 0.349, 0.384],
  "blue1livecom": [0.153, 0.231, 0.606],
  "blue2livecom": [0.033, 0.7, 0.97],
  "graylivecom": [0.427, 0.431, 0.439],
  "neoncyan": [0.0 , 0.9961, 0.988],
  "skobeloffgreen": [0, 0.455, 0.455],
  "tealblue": [0.211, 0.458, 0.533],
  "midnightgreen": [0, 0.286, 0.325],
  "turquoiseblue": [0, 1, 0.937],
  "teal": [0, 0.5, 0.5],
  "tealmodified": [0.2, 0.5, 0.5],
  "peacockblue": [0, 0.28, 0.34],
  "bluegreen": [0, 0.28, 0.34],
  "bluegreenpastel": [0.6, 0.87, 0.85],
  "bluepastel": [0.6, 0.6, 0.85],
  "lightturquoise": [0.63, 0.88, 0.88],
  "lightturquoise2": [0.58, 0.83, 0.83],
  "greenish": [0.54, 0.94, 0.89],
  "unisonblue": [0.58, 0.89, 0.96],
  "unisoncyan": [0.58, 0.89, 0.89],
}

def mygradient(N, color1, color2, final_value=False):
    """Generate a color gradient from color1 to color2"""
    if final_value:
        R = np.linspace(color1[0], color2[0], N)
        G = np.linspace(color1[1], color2[1], N)
        B = np.linspace(color1[2], color2[2], N)
    else:
        R = np.linspace(color1[0], color2[0], N)[1:-1]
        G = np.linspace(color1[1], color2[1], N)[1:-1]
        B = np.linspace(color1[2], color2[2], N)[1:-1]
    return np.vstack([R, G, B]).T
