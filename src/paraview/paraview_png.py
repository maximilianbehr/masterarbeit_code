# ### import the simple module from the paraview
from src.paraview.paraview.simple import *
# ### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
velopvd = PVDReader(
    FileName='/Users/daniels/Documents/LiClipseWorkspace/master/src/master/results/karman/chorin/1.3.0+/ref_1/velo.pvd')

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# set active source
SetActiveSource(velopvd)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [634, 297]

# show data in view
velopvdDisplay = Show(velopvd, renderView1)
# trace defaults for the display properties.
velopvdDisplay.ColorArrayName = [None, '']
velopvdDisplay.ScalarOpacityUnitDistance = 0.5614693652536232

# reset view to fit data
renderView1.ResetCamera()

# set scalar coloring
ColorBy(velopvdDisplay, ('POINTS', 'f_7'))

# rescale color and/or opacity maps used to include current data range
velopvdDisplay.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
velopvdDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'f7'
f7LUT = GetColorTransferFunction('f7')

# get opacity transfer function/opacity map for 'f7'
f7PWF = GetOpacityTransferFunction('f7')

# current camera placement for renderView1
renderView1.CameraPosition = [2.0, 0.5, 7.96522841660369]
renderView1.CameraFocalPoint = [2.0, 0.5, 0.0]
renderView1.CameraParallelScale = 2.0615528128088303

# current camera placement for renderView1
renderView1.CameraPosition = [2.0, 0.5, 7.96522841660369]
renderView1.CameraFocalPoint = [2.0, 0.5, 0.0]
renderView1.CameraParallelScale = 2.0615528128088303

# save animation images/movie
WriteAnimation('/Users/daniels/Desktop/vpng.png', Magnification=1, FrameRate=15.0, Compression=True)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [2.0, 0.5, 7.96522841660369]
renderView1.CameraFocalPoint = [2.0, 0.5, 0.0]
renderView1.CameraParallelScale = 2.0615528128088303

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).