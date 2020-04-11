from json import loads as json_loads
from json import dumps as json_dumps
from redis import Redis

class RedisDict:
    """A redis based dict."""

    def __init__(self, name, redis):
        self.name = name
        self.redis = redis

    def dict(self):
        return json_loads(redis.hget(self.name))

    def keys(self):
        return self.dict().keys()

    def __getitem__(self, key):
        if self.redis.exists(self.name):
            try:
                return json_loads(self.redis.hget(self.name, key))
            except TypeError:
                pass
            return None

    def __setitem__(self, key, value):
        self.redis.hset(self.name, key, json_dumps(value))

    def __delitem__(self, key):
        self.redis.hdel(self.name, key)
