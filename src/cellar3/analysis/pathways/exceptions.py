from src.cellar3.exceptions import ExternalResourceException


class PathwayException(ExternalResourceException):
    pass


class PathwayUnavailableException(PathwayException):
    pass


class ParsingException(PathwayException):
    pass


class GeneException(ParsingException):
    pass
