import originpro as op
import numpy as np
import sys
import matplotlib.pyplot as plt
from pathlib import Path

def test1():
    # Very useful, especially during development, when you are
    # liable to have a few uncaught exceptions.
    # Ensures that the Origin instance gets shut down properly.
    # Note: only applicable to external Python.

    def origin_shutdown_exception_hook(exctype, value, traceback):
        '''Ensures Origin gets shut down if an uncaught exception'''
        op.exit()
        sys.__excepthook__(exctype, value, traceback)
    if op and op.oext:
        sys.excepthook = origin_shutdown_exception_hook


    # Set Origin instance visibility.
    # Important for only external Python.
    # Should not be used with embedded Python. 
    if op.oext:
        op.set_show(True)


    # Example of plotting into a graph template.

    # Generate some data.
    arr_x = np.arange(100)
    arr_y1 = np.random.exponential(scale=2, size=100)
    arr_y2 = np.random.exponential(scale=3, size=100)

    # Create a new book with sheet and add generated data.
    wks = op.new_sheet()
    wks.from_list(0, arr_x, lname='X Values')
    wks.from_list(1, arr_y1, lname='Y1 Values')
    wks.from_list(2, arr_y2, lname='Y2 Values')

    # Create a new graph object based on a graph template.
    # Provide path to the template. We'll use a built-in
    # one but you can certainly use your own!
    gr = op.new_graph(template=op.path('e') + 'doubley.otp')

    # Add plots with data from worksheet to layer 1 and layer 2
    # and then rescale.
    gl_1 = gr[0]
    p1 = gl_1.add_plot(wks, coly='B', colx='A', type=202) # X is col A, Y is col B. 202 is Line + Symbol.
    p1.color = '#335eff'
    gl_1.rescale()

    gl_2 = gr[1]
    p2 = gl_2.add_plot(wks, 'C', 'A', type=202) # X is col A, Y is col C. 202 is Line + Symbol.
    p2.color = '#ff5833'
    gl_2.rescale()

    # change font size of axis titles
    gl_1.label('x', '@john:X Values')

    # Save the opju to your UFF.
    op.save(r'c:\tmp\test3.opju')


    # Exit running instance of Origin.
    # Required for external Python but don't use with embedded Python.
    if op.oext:
        op.exit()
        
def test2():


    # Create your plot
    plot = op.new_plot()

    # Access the X and Y axes
    x_axis = plot.xaxis
    y_axis = plot.yaxis

    # Set the desired font properties
    x_axis.font = "Arial"  # Change "Arial" to your preferred font
    x_axis.fontsize = 12

    y_axis.font = "Times New Roman"  # Change "Times New Roman" to your preferred font
    y_axis.fontsize = 14

    # Add your data and other plot elements...

    # Plot the graph
    op.plot(plot)
    
    
########################

def test3(filename):
    import win32com.client

    # Replace 'path\to\your\file.opju' with the actual path to your OPJU file
    opju_file = filename

    # Create Origin instance (replace 'Origin.Application' if necessary)
    origin = win32com.client.Dispatch('Origin.Application')

    # Open the OPJU file
    origin.FileOpen(opju_file)

    # Access the active worksheet (assuming the plot is in the first sheet)
    ws = origin.ActiveWorksheet

    # Get the active layer (assuming there's only one layer)
    layer = ws.Layers[0]

    # Access the X-axis object
    x_axis = layer.XAxis

    # Set the desired font properties (replace with your preferences)
    x_axis.Font = "Arial"  # Change "Arial" to your desired font
    x_axis.FontSize = 12

    # Close the workbook (optional)
    origin.Workbooks[0].Close()

    # Quit Origin (optional)
    origin.Quit()
    
def test4(filename):
    filename = Path(filename).with_suffix('.svg')
    x = np.linspace(0, 10, 100)
    y = np.sin(x)

    plt.scatter(x, y)
    plt.savefig(filename)
    print(f"Saved plot to {filename}")
    return filename
    
# import xml.etree.ElementTree as ET

# def extract_datapoints(svg_file):
#   """
#   Extracts datapoints (x, y coordinates) from a Matplotlib-generated SVG file.

#   Args:
#       svg_file (str): Path to the SVG file.

#   Returns:
#       list: A list of tuples containing (x, y) coordinates.
#   """
#   datapoints = []
#   tree = ET.parse(svg_file)
#   root = tree.getroot()
  
#   print(tree)

#   # Assuming datapoints are stored in path elements with specific attributes
#   # (modify these conditions based on your specific Matplotlib plot)
#   for path in root.findall(".//{http://www.w3.org/2000/svg}path[@d]"):
#     d_attribute = path.attrib["d"]
#     # Extract datapoints from the 'd' attribute using regular expressions (adapt based on format)
#     import re
#     matches = re.findall(r"[-\d+\.]+[ ,]", d_attribute)
#     for i in range(0, len(matches), 2):
#       x, y = matches[i], matches[i+1]
#       datapoints.append((float(x), float(y)))
      
#     # for d in datapoints:
#     #     print(d)

#   return np.asarray(datapoints)


import matplotlib.path as mpath
import xml.etree.ElementTree as ET

def extract_datapoints(svg_file):
  """
  Extracts datapoints (x, y) coordinates from a Matplotlib-generated SVG file.

  Args:
      svg_file (str): Path to the SVG file.

  Returns:
      list: A list of tuples containing (x, y) coordinates.
  """
  datapoints = []
  tree = ET.parse(svg_file)
  root = tree.getroot()

  # Assuming datapoints are stored in path elements
  for path in root.findall(".//{http://www.w3.org/2000/svg}path[@d]"):
    d_attribute = path.attrib["d"]
    # path_data = mpath.Path(d_attribute)  # Parse path data using Matplotlib's path library
    
    points = d_attribute.split(' ')[1:-1]  # Remove the 'M' and 'z' and split the string into a list
    vertices = [(float(points[i]), float(points[i+1])) for i in range(0, len(points), 2)]  # Convert the list into a list of tuples
    path_data = mpath.Path(vertices)
    
    for segment in path_data.itercurves():
      # Extract vertices (points) from each curve segment
      for vertex in segment.vertices:
        datapoints.append(vertex)

  return datapoints

filename = r'c:\tmp\test4.opju'
filename = test4(filename)
data= extract_datapoints(filename)
# print(data)
plt.scatter(data[:,0], data[:,1])
plt.show()
    