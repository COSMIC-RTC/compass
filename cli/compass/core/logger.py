"""
Logging configuration for COMPASS.

Provides consistent logging across all modules with Rich console output.
"""

import logging
from rich.logging import RichHandler
from rich.console import Console


def setup_logger(
    name: str = "compass",
    level: int = logging.INFO,
    verbose: bool = False,
) -> logging.Logger:
    """
    Set up a logger with Rich console output.
    
    Args:
        name: Logger name
        level: Logging level
        verbose: Enable verbose (DEBUG) logging
    
    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(name)
    
    # Set level
    if verbose:
        level = logging.DEBUG
    logger.setLevel(level)
    
    # Remove existing handlers
    logger.handlers.clear()
    
    # Add Rich handler for console output
    console = Console(stderr=True)
    rich_handler = RichHandler(
        console=console,
        rich_tracebacks=True,
        tracebacks_show_locals=verbose,
        markup=True,
    )
    rich_handler.setLevel(level)
    
    # Format
    formatter = logging.Formatter(
        "%(message)s",
        datefmt="[%X]",
    )
    rich_handler.setFormatter(formatter)
    logger.addHandler(rich_handler)
    
    return logger


class LoggerMixin:
    """Mixin class to add logging capability to any class."""
    
    @property
    def logger(self) -> logging.Logger:
        """Get logger for this class."""
        if not hasattr(self, '_logger'):
            self._logger = logging.getLogger(
                f"{self.__class__.__module__}.{self.__class__.__name__}"
            )
        return self._logger
