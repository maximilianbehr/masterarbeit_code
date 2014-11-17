#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
velopvd = PVDReader(FileName='/Users/daniels/Documents/LiClipseWorkspace/master/src/master/results/karman/chorin/1.3.0+/ref_1/velo.pvd')

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
velopvdDisplay.Representation = 'Surface'
velopvdDisplay.Ambient = 0.0
velopvdDisplay.AmbientColor = [1.0, 1.0, 1.0]
velopvdDisplay.ColorArrayName = [None, '']
velopvdDisplay.Diffuse = 1.0
velopvdDisplay.DiffuseColor = [1.0, 1.0, 1.0]
velopvdDisplay.LookupTable = None
velopvdDisplay.MapScalars = 1
velopvdDisplay.InterpolateScalarsBeforeMapping = 1
velopvdDisplay.Opacity = 1.0
velopvdDisplay.PointSize = 2.0
velopvdDisplay.LineWidth = 1.0
velopvdDisplay.Interpolation = 'Gouraud'
velopvdDisplay.Specular = 0.0
velopvdDisplay.SpecularColor = [1.0, 1.0, 1.0]
velopvdDisplay.SpecularPower = 100.0
velopvdDisplay.EdgeColor = [0.0, 0.0, 0.5]
velopvdDisplay.BackfaceRepresentation = 'Follow Frontface'
velopvdDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
velopvdDisplay.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
velopvdDisplay.BackfaceOpacity = 1.0
velopvdDisplay.Position = [0.0, 0.0, 0.0]
velopvdDisplay.Scale = [1.0, 1.0, 1.0]
velopvdDisplay.Orientation = [0.0, 0.0, 0.0]
velopvdDisplay.Origin = [0.0, 0.0, 0.0]
velopvdDisplay.Pickable = 1
velopvdDisplay.Texture = None
velopvdDisplay.NonlinearSubdivisionLevel = 1
velopvdDisplay.CubeAxesColor = [1.0, 1.0, 1.0]
velopvdDisplay.CubeAxesCornerOffset = 0.0
velopvdDisplay.CubeAxesFlyMode = 'Closest Triad'
velopvdDisplay.CubeAxesInertia = 1
velopvdDisplay.CubeAxesTickLocation = 'Inside'
velopvdDisplay.CubeAxesXAxisMinorTickVisibility = 1
velopvdDisplay.CubeAxesXAxisTickVisibility = 1
velopvdDisplay.CubeAxesXAxisVisibility = 1
velopvdDisplay.CubeAxesXGridLines = 0
velopvdDisplay.CubeAxesXTitle = 'X-Axis'
velopvdDisplay.CubeAxesUseDefaultXTitle = 1
velopvdDisplay.CubeAxesYAxisMinorTickVisibility = 1
velopvdDisplay.CubeAxesYAxisTickVisibility = 1
velopvdDisplay.CubeAxesYAxisVisibility = 1
velopvdDisplay.CubeAxesYGridLines = 0
velopvdDisplay.CubeAxesYTitle = 'Y-Axis'
velopvdDisplay.CubeAxesUseDefaultYTitle = 1
velopvdDisplay.CubeAxesZAxisMinorTickVisibility = 1
velopvdDisplay.CubeAxesZAxisTickVisibility = 1
velopvdDisplay.CubeAxesZAxisVisibility = 1
velopvdDisplay.CubeAxesZGridLines = 0
velopvdDisplay.CubeAxesZTitle = 'Z-Axis'
velopvdDisplay.CubeAxesUseDefaultZTitle = 1
velopvdDisplay.CubeAxesGridLineLocation = 'All Faces'
velopvdDisplay.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
velopvdDisplay.CustomBoundsActive = [0, 0, 0]
velopvdDisplay.OriginalBoundsRangeActive = [0, 0, 0]
velopvdDisplay.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
velopvdDisplay.CustomRangeActive = [0, 0, 0]
velopvdDisplay.UseAxesOrigin = 0
velopvdDisplay.AxesOrigin = [0.0, 0.0, 0.0]
velopvdDisplay.CubeAxesXLabelFormat = '%-#6.3g'
velopvdDisplay.CubeAxesYLabelFormat = '%-#6.3g'
velopvdDisplay.CubeAxesZLabelFormat = '%-#6.3g'
velopvdDisplay.StickyAxes = 0
velopvdDisplay.CenterStickyAxes = 0
velopvdDisplay.SelectionCellLabelBold = 0
velopvdDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
velopvdDisplay.SelectionCellLabelFontFamily = 'Arial'
velopvdDisplay.SelectionCellLabelFontSize = 18
velopvdDisplay.SelectionCellLabelItalic = 0
velopvdDisplay.SelectionCellLabelJustification = 'Left'
velopvdDisplay.SelectionCellLabelOpacity = 1.0
velopvdDisplay.SelectionCellLabelShadow = 0
velopvdDisplay.SelectionPointLabelBold = 0
velopvdDisplay.SelectionPointLabelColor = [0.5, 0.5, 0.5]
velopvdDisplay.SelectionPointLabelFontFamily = 'Arial'
velopvdDisplay.SelectionPointLabelFontSize = 18
velopvdDisplay.SelectionPointLabelItalic = 0
velopvdDisplay.SelectionPointLabelJustification = 'Left'
velopvdDisplay.SelectionPointLabelOpacity = 1.0
velopvdDisplay.SelectionPointLabelShadow = 0
velopvdDisplay.ScalarOpacityUnitDistance = 0.5614693652536232
velopvdDisplay.SelectMapper = 'Projected tetra'

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
f7LUT.InterpretValuesAsCategories = 0
f7LUT.EnableOpacityMapping = 0
f7LUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.0010825603718427003, 0.865003, 0.865003, 0.865003, 0.0021651207436854006, 0.705882, 0.0156863, 0.14902]
f7LUT.UseLogScale = 0
f7LUT.LockScalarRange = 0
f7LUT.ColorSpace = 'Diverging'
f7LUT.NanColor = [0.5, 0.0, 0.0]
f7LUT.Discretize = 1
f7LUT.NumberOfTableValues = 256
f7LUT.ScalarRangeInitialized = 1.0
f7LUT.HSVWrap = 0
f7LUT.VectorComponent = 0
f7LUT.VectorMode = 'Magnitude'
f7LUT.AllowDuplicateScalars = 1
f7LUT.Annotations = []
f7LUT.IndexedColors = []

# get opacity transfer function/opacity map for 'f7'
f7PWF = GetOpacityTransferFunction('f7')
f7PWF.Points = [0.0, 0.0, 0.5, 0.0, 0.0021651207436854006, 1.0, 0.5, 0.0]
f7PWF.AllowDuplicateScalars = 1
f7PWF.ScalarRangeInitialized = 1

# Rescale transfer function
f7LUT.RescaleTransferFunction(0.0, 220.388964926)

# Rescale transfer function
f7PWF.RescaleTransferFunction(0.0, 220.388964926)

# Properties modified on f7LUT
f7LUT.LockScalarRange = 1

# Properties modified on animationScene1
animationScene1.PlayMode = 'Real Time'

# Properties modified on animationScene1
animationScene1.Duration = 8

# current camera placement for renderView1
renderView1.CameraPosition = [2.0, 0.5, 7.96522841660369]
renderView1.CameraFocalPoint = [2.0, 0.5, 0.0]
renderView1.CameraParallelScale = 2.0615528128088303

# current camera placement for renderView1
renderView1.CameraPosition = [2.0, 0.5, 7.96522841660369]
renderView1.CameraFocalPoint = [2.0, 0.5, 0.0]
renderView1.CameraParallelScale = 2.0615528128088303

# save animation images/movie
WriteAnimation('/Applications/v.avi', Magnification=1, FrameRate=15.0, Compression=True)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [2.0, 0.5, 7.96522841660369]
renderView1.CameraFocalPoint = [2.0, 0.5, 0.0]
renderView1.CameraParallelScale = 2.0615528128088303

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).