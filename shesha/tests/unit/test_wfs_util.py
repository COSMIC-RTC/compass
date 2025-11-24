"""
Tests for shesha.util.wfs_util module
"""

import numpy as np
import pytest
import tempfile
import os
from shesha.util.wfs_util import write_wfs_custom_fits, add_doc_content


class TestWriteWFSCustomFits:
    """Test WFS FITS file writing."""
    
    def test_write_wfs_custom_fits_basic(self):
        """Test basic WFS FITS file writing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "test_wfs.fits")
            
            xpos = np.array([0.0, 1.0, 2.0], dtype=np.float32)
            ypos = np.array([0.0, 1.0, 2.0], dtype=np.float32)
            
            hdul = write_wfs_custom_fits(
                filepath,
                "sh",
                xpos,
                ypos,
                xcenter=128,
                ycenter=128,
                pixsize=0.01,
                pupm=8.0,
                subap_diam=0.5,
            )
            
            assert os.path.exists(filepath)
            assert hdul is not None
            assert len(hdul) == 2
    
    def test_write_wfs_custom_fits_header_content(self):
        """Test WFS FITS header contains correct information."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "test_wfs.fits")
            
            xpos = np.array([0.0, 1.0], dtype=np.float32)
            ypos = np.array([0.0, 1.0], dtype=np.float32)
            xcenter = 100.5
            ycenter = 101.5
            pixsize = 0.012
            pupm = 10.0
            subap_diam = 0.6
            
            hdul = write_wfs_custom_fits(
                filepath,
                "sh",
                xpos,
                ypos,
                xcenter=xcenter,
                ycenter=ycenter,
                pixsize=pixsize,
                pupm=pupm,
                subap_diam=subap_diam,
            )
            
            header = hdul[0].header
            
            assert header["TYPE"] == "sh"
            assert header["XCENTER"] == xcenter
            assert header["YCENTER"] == ycenter
            assert np.isclose(header["PIXSIZE"], pixsize)
            assert np.isclose(header["PUPM"], pupm)
            assert np.isclose(header["SUBAPD"], subap_diam)
    
    def test_write_wfs_custom_fits_data_extension(self):
        """Test WFS FITS data extension."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "test_wfs.fits")
            
            xpos = np.array([0.0, 1.0, 2.0], dtype=np.float32)
            ypos = np.array([0.5, 1.5, 2.5], dtype=np.float32)
            
            hdul = write_wfs_custom_fits(
                filepath,
                "sh",
                xpos,
                ypos,
                xcenter=128,
                ycenter=128,
                pixsize=0.01,
                pupm=8.0,
                subap_diam=0.5,
            )
            
            # Check image extension
            image_hdu = hdul[1]
            assert image_hdu.name == "XPOS_YPOS"
            assert image_hdu.data.shape == (2, 3)
    
    def test_write_wfs_custom_fits_invalid_type(self):
        """Test WFS FITS with invalid WFS type."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "test_wfs.fits")
            
            xpos = np.array([0.0], dtype=np.float32)
            ypos = np.array([0.0], dtype=np.float32)
            
            with pytest.raises(RuntimeError):
                write_wfs_custom_fits(
                    filepath,
                    "invalid_type",
                    xpos,
                    ypos,
                    xcenter=128,
                    ycenter=128,
                    pixsize=0.01,
                    pupm=8.0,
                    subap_diam=0.5,
                )
    
    def test_write_wfs_custom_fits_large_array(self):
        """Test WFS FITS with large subaperture array."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "test_wfs.fits")
            
            n_subap = 100
            xpos = np.linspace(0, 10, n_subap, dtype=np.float32)
            ypos = np.linspace(0, 10, n_subap, dtype=np.float32)
            
            hdul = write_wfs_custom_fits(
                filepath,
                "sh",
                xpos,
                ypos,
                xcenter=128,
                ycenter=128,
                pixsize=0.01,
                pupm=8.0,
                subap_diam=0.5,
            )
            
            assert hdul[1].data.shape == (2, n_subap)
    
    def test_write_wfs_custom_fits_small_array(self):
        """Test WFS FITS with single subaperture."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "test_wfs.fits")
            
            xpos = np.array([0.0], dtype=np.float32)
            ypos = np.array([0.0], dtype=np.float32)
            
            hdul = write_wfs_custom_fits(
                filepath,
                "sh",
                xpos,
                ypos,
                xcenter=128,
                ycenter=128,
                pixsize=0.01,
                pupm=8.0,
                subap_diam=0.5,
            )
            
            assert hdul[1].data.shape == (2, 1)
    
    def test_write_wfs_custom_fits_preserves_data(self):
        """Test that written data matches input data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "test_wfs.fits")
            
            xpos = np.array([0.1, 0.5, 1.2, 2.3], dtype=np.float32)
            ypos = np.array([0.2, 0.6, 1.3, 2.4], dtype=np.float32)
            
            hdul = write_wfs_custom_fits(
                filepath,
                "sh",
                xpos,
                ypos,
                xcenter=128,
                ycenter=128,
                pixsize=0.01,
                pupm=8.0,
                subap_diam=0.5,
            )
            
            # Extract data from FITS file (data stored as [xpos, ypos])
            data = hdul[1].data
            xpos_read = data[0, :]
            ypos_read = data[1, :]
            
            assert np.allclose(xpos_read, xpos, rtol=1e-5)
            assert np.allclose(ypos_read, ypos, rtol=1e-5)
    
    def test_write_wfs_custom_fits_coordinates_as_float64(self):
        """Test that coordinates are stored as float64 in FITS."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "test_wfs.fits")
            
            xpos = np.array([0.0, 1.0], dtype=np.float32)
            ypos = np.array([0.0, 1.0], dtype=np.float32)
            
            hdul = write_wfs_custom_fits(
                filepath,
                "sh",
                xpos,
                ypos,
                xcenter=128,
                ycenter=128,
                pixsize=0.01,
                pupm=8.0,
                subap_diam=0.5,
            )
            
            assert hdul[1].data.dtype == np.float64
    
    def test_write_wfs_custom_fits_overwrite(self):
        """Test that existing files are overwritten."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "test_wfs.fits")
            
            # Write first file
            xpos1 = np.array([0.0], dtype=np.float32)
            ypos1 = np.array([0.0], dtype=np.float32)
            write_wfs_custom_fits(
                filepath,
                "sh",
                xpos1,
                ypos1,
                xcenter=100,
                ycenter=100,
                pixsize=0.01,
                pupm=8.0,
                subap_diam=0.5,
            )
            
            # Write second file with different data
            xpos2 = np.array([1.0, 2.0, 3.0], dtype=np.float32)
            ypos2 = np.array([1.0, 2.0, 3.0], dtype=np.float32)
            hdul = write_wfs_custom_fits(
                filepath,
                "sh",
                xpos2,
                ypos2,
                xcenter=200,
                ycenter=200,
                pixsize=0.02,
                pupm=10.0,
                subap_diam=0.6,
            )
            
            # Check that second data is used
            assert hdul[0].header["XCENTER"] == 200
            assert hdul[1].data.shape == (2, 3)
    
    def test_write_wfs_custom_fits_negative_coordinates(self):
        """Test WFS FITS with negative coordinates."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "test_wfs.fits")
            
            xpos = np.array([-1.0, 0.0, 1.0], dtype=np.float32)
            ypos = np.array([-2.0, 0.0, 2.0], dtype=np.float32)
            
            hdul = write_wfs_custom_fits(
                filepath,
                "sh",
                xpos,
                ypos,
                xcenter=128,
                ycenter=128,
                pixsize=0.01,
                pupm=8.0,
                subap_diam=0.5,
            )
            
            # Data should be preserved with sign (data stored as [xpos, ypos])
            data = hdul[1].data
            assert np.isclose(data[0, 0], -1.0)
            assert np.isclose(data[1, 0], -2.0)


class TestAddDocContent:
    """Test documentation decorator."""
    
    def test_add_doc_content_decorator(self):
        """Test that add_doc_content decorator works."""
        test_content = "test documentation"
        
        @add_doc_content(test_content)
        def dummy_func():
            """Function with {0}"""
            pass
        
        assert test_content in dummy_func.__doc__
    
    def test_add_doc_content_preserves_function(self):
        """Test that decorator preserves function behavior."""
        @add_doc_content("doc")
        def test_func(x):
            """Function with {0}"""
            return x * 2
        
        assert test_func(5) == 10


class TestWFSIntegration:
    """Integration tests for WFS utilities."""
    
    def test_wfs_fits_roundtrip(self):
        """Test writing and reading WFS FITS file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "test_wfs.fits")
            
            xpos_orig = np.array([0.1, 0.5, 1.0, 2.5], dtype=np.float32)
            ypos_orig = np.array([0.2, 0.6, 1.1, 2.6], dtype=np.float32)
            
            write_wfs_custom_fits(
                filepath,
                "sh",
                xpos_orig,
                ypos_orig,
                xcenter=128.5,
                ycenter=129.5,
                pixsize=0.0125,
                pupm=8.5,
                subap_diam=0.55,
            )
            
            # Read back the file
            from astropy.io import fits
            with fits.open(filepath) as hdul:
                data = hdul[1].data
                xpos_read = data[0, :]
                ypos_read = data[1, :]
                
                assert np.allclose(xpos_read, xpos_orig, rtol=1e-5)
                assert np.allclose(ypos_read, ypos_orig, rtol=1e-5)
    
    def test_wfs_fits_multiple_files(self):
        """Test creating multiple WFS FITS files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            files = []
            
            for i in range(3):
                filepath = os.path.join(tmpdir, f"wfs_{i}.fits")
                xpos = np.linspace(0, i, i+2, dtype=np.float32)
                ypos = np.linspace(0, i, i+2, dtype=np.float32)
                
                write_wfs_custom_fits(
                    filepath,
                    "sh",
                    xpos,
                    ypos,
                    xcenter=100 + i,
                    ycenter=100 + i,
                    pixsize=0.01 + i*0.001,
                    pupm=8.0 + i,
                    subap_diam=0.5,
                )
                
                files.append(filepath)
            
            # Check all files exist
            for f in files:
                assert os.path.exists(f)
