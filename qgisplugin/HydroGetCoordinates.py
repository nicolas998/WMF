from qgis.gui import QgsMapTool

class PointTool(QgsMapTool):
    def __init__(self, canvas, spinLat, spinLon):
        QgsMapTool.__init__(self, canvas)
        self.canvas = canvas
        self.spinLat = spinLat
        self.spinLon = spinLon

    def canvasPressEvent(self, event):
        pass

    def canvasMoveEvent(self, event):
        x = event.pos().x()
        y = event.pos().y()
        #point = self.canvas.getCoordinateTransform().toMapCoordinates(x, y)
        #self.spinLat.setValue(point.y())
        #self.spinLon.setValue(point.x())

    def canvasReleaseEvent(self, event):
        #Get the click
        x = event.pos().x()
        y = event.pos().y()
        point = self.canvas.getCoordinateTransform().toMapCoordinates(x, y)
        self.spinLat.setValue(point.y())
        self.spinLon.setValue(point.x())

    def activate(self):
        pass

    def deactivate(self):
        pass

    def isZoomTool(self):
        return False

    def isTransient(self):
        return False

    def isEditTool(self):
        return True
