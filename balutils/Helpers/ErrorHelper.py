# Central helper class for error messages

class ErrorHelper:

    @classmethod
    def missing_column(cls, col):
        return AttributeError(
            "{} not found in joined ".format(col) + "catalog but required for requested cuts!"
        )