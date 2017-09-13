from datetime import timedelta, tzinfo


class FixedOffset(tzinfo):
    """Fixed offset in minutes: `time = utc_time + utc_offset`.

    Inspriation: https://stackoverflow.com/questions/1101508/how-to-parse-dates-with-0400-timezone-string-in-python"""
    def __init__(self, offset):
        self.__offset = timedelta(minutes=offset)
        hours, minutes = divmod(offset, 60)
        self.__name = '<%+03d%02d>%+d' % (hours, minutes, -hours)

    def utcoffset(self, dt=None):
        return self.__offset

    def tzname(self, dt=None):
        return self.__name

    def dst(self, dt=None):
        return timedelta(0)

    def __repr__(self):
        return 'FixedOffset(%d)' % (self.utcoffset().total_seconds() / 60)
