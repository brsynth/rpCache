from json import loads as json_loads
from json import dumps as json_dumps
from redis import Redis


class RedisDict:
    """A redis based dict."""

    def __init__(self, name, redis, data={}):
        self.name = name
        self.redis = redis
        # This avoids to load all dict from redis each time we access to a key (often). So better than 'hmset'
        for key in data:
            self.__setitem__(key, data[key])

    def dict(self):
        return self.redis.hgetall(self.name)

    def keys(self):
        return self.redis.hkeys(self.name)

    def exists(self):
        return self.redis.exists(self.name)

    def len(self):
        return self.redis.hlen(self.name)

    def is_empty(self):
        return self.len()==0

    def __iter__(self):
        return iter(self.keys())

    def __contains__(self, key):
        return self.redis.hexists(self.name, key)

    def __getitem__(self, key):
        item = self.redis.hget(self.name, key)
        # JSON for nested dictionnaries
        if item: return json_loads(item)
        else: raise KeyError

    def __setitem__(self, key, value):
        # JSON for nested dictionnaries
        self.redis.hset(self.name, key, json_dumps(value))

    def __eq__(self, redis_dict):
        self.name = redis_dict.name
        self.redis = redis_dict.redis
        for key in redis_dict.keys():
            self.__setitem__(key, redis_dict[key])

    def update(self, redis_dict):
        for field in redis_dict:
            self.__setitem__(field, redis_dict[field])
