#!/usr/env/python

## Import General Tools
import sys
import os
import argparse
import logging

from ginga import AstroImage, colors
from ginga.misc import log
from ginga.qtw.QtHelp import QtGui, QtCore
from ginga.qtw.ImageViewQt import CanvasView, ScrolledView
from ginga.canvas.CanvasObject import get_canvas_types

from pathlib import Path
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord, Angle
from astropy.wcs import WCS, FITSFixedWarning
from astropy.visualization import MinMaxInterval, PercentileInterval, ImageNormalize
from astroquery.vizier import Vizier

from scipy.signal import medfilt

Vizier.ROW_LIMIT = 1000


##-------------------------------------------------------------------------
## Parse Command Line Arguments
##-------------------------------------------------------------------------
## create a parser object for understanding command-line arguments
p = argparse.ArgumentParser(description='''
''')
## add flags
p.add_argument("-v", "--verbose", dest="verbose",
    default=False, action="store_true",
    help="Be verbose! (default = False)")
## add options
p.add_argument("--catalog", "-c", dest="catalog", type=str,
    default="UCAC5",
    help="The name of the stellar catalog to use.")
## add arguments
p.add_argument('files', nargs='*',
               help="Files to check.")
args = p.parse_args()


##-------------------------------------------------------------------------
## Create logger object
##-------------------------------------------------------------------------
log = logging.getLogger('CheckWCS')
log.setLevel(logging.INFO)
if args.verbose is True:
    log.setLevel(logging.DEBUG)
## Set up console output
LogConsoleHandler = logging.StreamHandler()
LogConsoleHandler.setLevel(logging.DEBUG)
LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s',
                              datefmt='%Y-%m-%d %H:%M:%S')
LogConsoleHandler.setFormatter(LogFormat)
log.addHandler(LogConsoleHandler)
## Set up file output
# LogFileName = None
# LogFileHandler = logging.FileHandler(LogFileName)
# LogFileHandler.setLevel(logging.DEBUG)
# LogFileHandler.setFormatter(LogFormat)
# log.addHandler(LogFileHandler)


##-------------------------------------------------------------------------
## Ginga FITS Viewer
##-------------------------------------------------------------------------
class FitsViewer(QtGui.QMainWindow):

    def __init__(self, logger):
        super(FitsViewer, self).__init__()
        self.logger = logger
        self.drawcolors = colors.get_colors()
        self.dc = get_canvas_types()

        self.set_medfilt(0)

        # create the ginga viewer and configure it
        fi = CanvasView(self.logger, render='widget')
        fi.enable_autocuts('on')
        fi.set_autocut_params('zscale')
        fi.enable_autozoom('on')
        fi.set_callback('drag-drop', self.drop_file)
        fi.add_callback('cursor-changed', self.cursor_cb)
        fi.set_bg(0.2, 0.2, 0.2)
        fi.ui_set_active(True)
        self.fitsimage = fi

        # enable some user interaction
        bd = fi.get_bindings()
        bd.enable_all(True)

        # canvas that we will draw on
        canvas = self.dc.DrawingCanvas()
#         canvas.enable_draw(True)
#         canvas.enable_edit(True)
#         canvas.set_drawtype('circle', color='lightblue')
        canvas.set_surface(fi)
        self.canvas = canvas
        # add canvas to view
        #fi.add(canvas)
        private_canvas = fi.get_canvas()
        private_canvas.add(canvas)
#         canvas.register_for_cursor_drawing(fi)
#         canvas.add_callback('draw-event', self.draw_cb)
#         canvas.set_draw_mode('draw')
        canvas.ui_set_active(True)

        self.drawtypes = canvas.get_drawtypes()
        self.drawtypes.sort()

        w = fi.get_widget()
        w.resize(1000, 1000)

        # add scrollbar interface around this viewer
        si = ScrolledView(fi)

        vbox = QtGui.QVBoxLayout()
        vbox.setContentsMargins(QtCore.QMargins(2, 2, 2, 2))
        vbox.setSpacing(1)
        vbox.addWidget(si, stretch=1)

        self.readout = QtGui.QLabel("")
        vbox.addWidget(self.readout, stretch=0,
                       alignment=QtCore.Qt.AlignCenter)

        hbox = QtGui.QHBoxLayout()
        hbox.setContentsMargins(QtCore.QMargins(4, 2, 4, 2))

        wsolve = QtGui.QPushButton("Solve Astrometry")
        wsolve.clicked.connect(self.solve_astrometry)
        woverlay = QtGui.QPushButton("Overlay Catalog Stars")
        woverlay.clicked.connect(self.overlay_catalog)
        wclear = QtGui.QPushButton("Clear Overlays")
        wclear.clicked.connect(self.clear_overlays)

        wmedfilt = QtGui.QCheckBox("Median Filter")
        wmedfilt.stateChanged.connect(self.set_medfilt)
        self.wmedfilt = wmedfilt

        wopen = QtGui.QPushButton("Open File")
        wopen.clicked.connect(self.open_file)
        wquit = QtGui.QPushButton("Quit")
        wquit.clicked.connect(self.quit)

        hbox.addStretch(1)
        hbox.addWidget(wsolve, stretch=0)
        hbox.addWidget(woverlay, stretch=0)
        hbox.addWidget(wclear, stretch=0)
        hbox.addWidget(wmedfilt, stretch=0)
        hbox.addWidget(wopen, stretch=0)
        hbox.addWidget(wquit, stretch=0)

        hw = QtGui.QWidget()
        hw.setLayout(hbox)
        vbox.addWidget(hw, stretch=0)

        vw = QtGui.QWidget()
        self.setCentralWidget(vw)
        vw.setLayout(vbox)

    def set_medfilt(self, kind):
        self.medfilt = {0: False, 2: True}[kind]
        self.medfilt_str = {False: '', True: '(filtered)'}[self.medfilt]
        log.info(f"Setting median filter flag to {self.medfilt}")

    def cursor_cb(self, viewer, button, data_x, data_y):
        """This gets called when the data position relative to the cursor
        changes.
        """
        # Get the value under the data coordinates
        try:
            # We report the value across the pixel, even though the coords
            # change halfway across the pixel
            value = viewer.get_data(int(data_x + viewer.data_off),
                                    int(data_y + viewer.data_off))

        except Exception:
            value = None

        fits_x, fits_y = data_x + 1, data_y + 1

        # Calculate WCS RA
        try:
            # NOTE: image function operates on DATA space coords
            image = viewer.get_image()
            if image is None:
                # No image loaded
                return
            ra_txt, dec_txt = image.pixtoradec(fits_x, fits_y,
                                               format='str', coords='fits')
        except Exception as e:
            self.logger.warning("Bad coordinate conversion: %s" % (
                str(e)))
            ra_txt = 'BAD WCS'
            dec_txt = 'BAD WCS'

        text = "RA: %s  DEC: %s  X: %.2f  Y: %.2f  Value: %s" % (
            ra_txt, dec_txt, fits_x, fits_y, value)
        self.readout.setText(text)

    def load_file(self, filepath):
        self.clear_overlays()
        image = AstroImage.AstroImage(logger=self.logger)
        image.load_file(filepath)
        if self.medfilt is True:
            data = image.get_data()
            filtered_image = medfilt(data, (3,3))
            image.set_data(filtered_image)
        self.image = image
        self.fitsimage.set_image(image)
        self.setWindowTitle(f'{filepath} {self.medfilt_str}')
        self.overlay_catalog()

    def open_file(self):
        res = QtGui.QFileDialog.getOpenFileName(self, "Open FITS file",
                                                ".", "FITS files (*.fits)")
        if isinstance(res, tuple):
            fileName = res[0]
        else:
            fileName = str(res)
        if len(fileName) != 0:
            self.load_file(fileName)

    def drop_file(self, fitsimage, paths):
        fileName = paths[0]
        self.load_file(fileName)

    def quit(self, *args):
        self.logger.info("Attempting to shut down the application...")
        self.deleteLater()

    def clear_overlays(self):
        self.canvas.delete_all_objects()

    def solve_astrometry(self):
        log.warning('Astrometry solver not yet implemented.')

    def overlay_catalog(self):
        # Get WCS Info From Image Header
        h = self.image.get_header()
        wcs = WCS(h)

        nx, ny = self.image.shape
        result = wcs.all_pix2world(np.array(nx/2), np.array(ny/2), 1)
        center = SkyCoord(result[0], result[1], frame='icrs', unit='deg')
        center_str = center.to_string("hmsdms", precision=1, sep=':')
        footprint = wcs.calc_footprint()
        FoV = ((max(footprint[:,0])-min(footprint[:,0]))*np.cos(center.dec.deg*np.pi/180),
                max(footprint[:,1])-min(footprint[:,1]))
        log.info(f'Found WCS with center coordinate: {center_str}')
        log.info(f'WCS has approximate diameter of {max(FoV):.2f} deg '\
                 f'({max(FoV)*60:.1f} arcmin)')
        FoV = ((max(footprint[:,0])-min(footprint[:,0])),
                max(footprint[:,1])-min(footprint[:,1]))

        catalog_names = ['UCAC5', 'Gaia DR2', '2MASS PSC']
        catalog_IDs = ['I/340', 'I/345', 'II/246']
        colors = ['lightblue', 'green', 'red']
        radii = [10, 12, 14]
        RAcolname = ['RAJ2000', 'RA_ICRS', 'RAJ2000']
        DEcolname = ['DEJ2000', 'DE_ICRS', 'DEJ2000']

        result = Vizier.query_region(center,
                                     width=1.05*FoV[0]*u.deg,
                                     height=1.05*FoV[1]*u.deg,
                                     catalog=catalog_IDs)

        # Overlay Stars
        for i,name in enumerate(catalog_names):
            log.info(f'Retrieved {len(result[i])} stars from {name} catalog')
            for j,star in enumerate(result[i]):
                c = SkyCoord(star[RAcolname[i]]*u.deg, star[DEcolname[i]]*u.deg)
                c_str = c.to_string('hmsdms', sep=':', precision=1)
                x, y = wcs.all_world2pix(c.ra.deg, c.dec.deg, 1)
                log.debug(f"{name} {c_str} {x:.1f} {y:.1f}")
                if (x > 0) and (x <= nx) and (y > 0) and (y <= ny):
                    log.debug(f"Plotting {name}: {c_str} {x:.1f} {y:.1f}")
                    self.canvas.add(self.dc.Circle(x, y, radius=radii[i],
                                    alpha=0.7, color=colors[i]))



##-------------------------------------------------------------------------
## Main Program
##-------------------------------------------------------------------------
def main(file, ext=0):

    app = QtGui.QApplication(sys.argv)

    # Display File in Viewer Window
    w = FitsViewer(log)
    w.resize(1000, 1000)
    w.show()
    app.setActiveWindow(w)
    w.raise_()
    w.activateWindow()
    w.load_file(str(file))

    app.exec_()

if __name__ == '__main__':
    for file in args.files:
        f = Path(file).expanduser()
        main(f)
