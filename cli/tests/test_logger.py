"""
Tests for compass.core.logger module.
"""

import logging
from compass.core.logger import setup_logger, LoggerMixin


class TestSetupLogger:
    """Test cases for setup_logger function."""
    
    def test_setup_logger_default_params(self):
        """Test setup_logger with default parameters."""
        logger = setup_logger()
        
        assert isinstance(logger, logging.Logger)
        assert logger.name == "compass"
        assert logger.level == logging.INFO
    
    def test_setup_logger_custom_name(self):
        """Test setup_logger with custom name."""
        logger = setup_logger(name="custom_logger")
        
        assert logger.name == "custom_logger"
    
    def test_setup_logger_debug_level(self):
        """Test setup_logger with DEBUG level."""
        logger = setup_logger(level=logging.DEBUG)
        
        assert logger.level == logging.DEBUG
    
    def test_setup_logger_verbose_true(self):
        """Test setup_logger with verbose=True."""
        logger = setup_logger(verbose=True)
        
        # Verbose should set DEBUG level
        assert logger.level == logging.DEBUG
    
    def test_setup_logger_verbose_overrides_level(self):
        """Test setup_logger verbose overrides level parameter."""
        logger = setup_logger(level=logging.INFO, verbose=True)
        
        # Verbose should override and set DEBUG
        assert logger.level == logging.DEBUG
    
    def test_setup_logger_has_handlers(self):
        """Test setup_logger adds handlers."""
        logger = setup_logger()
        
        assert len(logger.handlers) > 0
    
    def test_setup_logger_multiple_calls_same_name(self):
        """Test setup_logger replaces handlers on multiple calls."""
        logger1 = setup_logger(name="test_logger")
        handler_count_1 = len(logger1.handlers)
        
        logger2 = setup_logger(name="test_logger")
        handler_count_2 = len(logger2.handlers)
        
        # Should clear handlers and add fresh ones
        assert logger1 is logger2
        assert handler_count_2 == handler_count_1


class TestLoggerMixin:
    """Test cases for LoggerMixin class."""
    
    def test_logger_mixin_has_logger_property(self):
        """Test LoggerMixin provides logger property."""
        
        class TestClass(LoggerMixin):
            pass
        
        obj = TestClass()
        assert hasattr(obj, 'logger')
        assert isinstance(obj.logger, logging.Logger)
    
    def test_logger_mixin_logger_name(self):
        """Test LoggerMixin logger has correct name."""
        
        class MyModule(LoggerMixin):
            pass
        
        obj = MyModule()
        logger = obj.logger
        
        # Logger name should be based on module and class
        assert "MyModule" in logger.name
    
    def test_logger_mixin_caches_logger(self):
        """Test LoggerMixin caches the logger instance."""
        
        class TestClass(LoggerMixin):
            pass
        
        obj = TestClass()
        logger1 = obj.logger
        logger2 = obj.logger
        
        assert logger1 is logger2
    
    def test_logger_mixin_different_instances_same_name(self):
        """Test different instances get same logger name."""
        
        class MyClass(LoggerMixin):
            pass
        
        obj1 = MyClass()
        obj2 = MyClass()
        
        assert obj1.logger.name == obj2.logger.name


class TestLoggerIntegration:
    """Integration tests for logging setup."""
    
    def test_setup_multiple_loggers(self):
        """Test setting up multiple loggers."""
        logger_a = setup_logger(name="logger_a")
        logger_b = setup_logger(name="logger_b")
        
        assert logger_a.name == "logger_a"
        assert logger_b.name == "logger_b"
        assert logger_a is not logger_b
    
    def test_mixin_with_setup_logger(self):
        """Test LoggerMixin works with setup_logger."""
        
        class MyService(LoggerMixin):
            def __init__(self):
                # Get logger for this service
                self._service_logger = setup_logger(
                    name=f"{self.__class__.__module__}.{self.__class__.__name__}"
                )
            
            @property
            def service_logger(self):
                return self._service_logger
        
        service = MyService()
        assert isinstance(service.service_logger, logging.Logger)
